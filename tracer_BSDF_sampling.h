// Last confirm date 2025/12/06 
#ifndef TRACER_BSDF_SAMPLING_H
#define TRACER_BSDF_SAMPLING_H

#include "Speed_up.h"
#include "tracer_light_source.h"

// ============================================================
// Phong specular 開關
//   1: 啟用 Phong 高光 (原本行為)
//   0: 關閉 Phong，高光能量併入 diffuse，整體亮度不變
// ============================================================
#ifndef ENABLE_PHONG_SPECULAR
#define ENABLE_PHONG_SPECULAR 0     // 1:on, 0:off
#endif

// Kd, Ks 的「有效值」：考慮 Phong switch 後的版本
inline void get_effective_kd_ks(const Material* mat,
                                double* kd_out,
                                double* ks_out) {
    double kd = mat->kd;
    double ks = mat->ks;
    if (kd < 0.0) kd = 0.0;
    if (ks < 0.0) ks = 0.0;

#if ENABLE_PHONG_SPECULAR
    // 正常情況：Kd, Ks 如材質設定
    *kd_out = kd;
    *ks_out = ks;
#else
    // 關掉 Phong：把 Ks 能量丟回 diffuse
    *kd_out = kd + ks;
    *ks_out = 0.0;
#endif
}

// ----------------- ONB -----------------

struct ONB {
    vec3 u, v, w;
};

inline ONB build_onb(const vec3& n) {
    ONB onb;
    onb.w = unit_vector(n);
    vec3 a;
    if (fabs(onb.w.x()) > 0.9) a = vec3(0,1,0);
    else a = vec3(1,0,0);
    onb.v = unit_vector(cross(onb.w, a));
    onb.u = cross(onb.v, onb.w);
    return onb;
}

inline vec3 local_to_world(const ONB& onb, const vec3& a) {
    return a.x() * onb.u + a.y() * onb.v + a.z() * onb.w;
}

// ----------------- BSDF sampling（半球 / Phong lobe） -----------------

inline vec3 sample_cosine_hemisphere(const vec3& n, double u1, double u2) {
    double r = sqrt(u1);
    double theta = 2.0 * M_PI * u2;
    double x = r * cos(theta);
    double y = r * sin(theta);
    double z = sqrt(fmax(0.0, 1.0 - u1));

    ONB onb = build_onb(n);
    return local_to_world(onb, vec3(x,y,z));
}

inline double cosine_hemisphere_pdf(const vec3& n, const vec3& wi) {
    vec3 wi_n = unit_vector(wi);
    double cos_theta = dot(n, wi_n);
    if (cos_theta < 0.0) cos_theta = 0.0;
    return cos_theta > 0.0 ? cos_theta / M_PI : 0.0;
}

inline vec3 reflect_vec(const vec3& v, const vec3& n) {
    return v - 2.0 * dot(v, n) * n;
}

inline vec3 sample_phong_lobe(const vec3& rdir, double shininess,
                              double u1, double u2) {
    double cos_alpha = pow(u1, 1.0 / (shininess + 1.0));
    double sin_alpha = sqrt(fmax(0.0, 1.0 - cos_alpha * cos_alpha));
    double phi = 2.0 * M_PI * u2;
    double x = sin_alpha * cos(phi);
    double y = sin_alpha * sin(phi);
    double z = cos_alpha;

    ONB onb = build_onb(rdir);
    return local_to_world(onb, vec3(x,y,z));
}

inline double phong_lobe_pdf(double shininess,
                             const vec3& rdir,
                             const vec3& wi) {
    vec3 r = unit_vector(rdir);
    vec3 wi_n = unit_vector(wi);
    double cos_alpha = dot(r, wi_n);
    if (cos_alpha < 0.0) cos_alpha = 0.0;
    if (cos_alpha <= 0.0) return 0.0;
    return (shininess + 1.0) * pow(cos_alpha, shininess) / (2.0 * M_PI);
}

// ----------------- BRDF eval（Lambert + Phong） -----------------

inline vec3 eval_brdf(const Material* mat,
                      const vec3& n,
                      const vec3& wi,
                      const vec3& wo) {
    vec3 wi_n = unit_vector(wi);
    vec3 wo_n = unit_vector(wo);

    double cos_theta_i = dot(n, wi_n);
    if (cos_theta_i < 0.0) cos_theta_i = 0.0;
    if (cos_theta_i <= 0.0) return vec3(0.0, 0.0, 0.0);

    double kd_eff, ks_eff;
    get_effective_kd_ks(mat, &kd_eff, &ks_eff);

    // Diffuse: ρ_d / π
    vec3 fd = (kd_eff / M_PI) * mat->color;

    // Specular (Phong)
    vec3 fs(0.0,0.0,0.0);
#if ENABLE_PHONG_SPECULAR
    if (ks_eff > 0.0) {
        vec3 rdir = reflect_vec(-wo_n, n);
        double cos_alpha = dot(unit_vector(rdir), wi_n);
        if (cos_alpha < 0.0) cos_alpha = 0.0;
        if (cos_alpha > 0.0) {
            double norm = (mat->shininess + 2.0) / (2.0 * M_PI);
            fs = ks_eff * norm * pow(cos_alpha, mat->shininess) * mat->color;
        }
    }
#endif

    return fd + fs;
}

// ----------------- BRDF direction sampling -----------------

inline int sample_brdf_direction(const Material* mat,
                                 const vec3& n,
                                 const vec3& wo,
                                 vec3* wi_out,
                                 double* pdf_out) {
    double kd, ks;
    get_effective_kd_ks(mat, &kd, &ks);

    double sum = kd + ks;
    double pd = 1.0, ps = 0.0;
    if (sum > 0.0) {
        pd = kd / sum;
        ps = ks / sum;
    }

    double u = rand01();
    int use_diffuse = (u < pd);

    vec3 wi;
    // ---- Diffuse branch ----
    if (use_diffuse || ks <= 0.0) {
        double u1 = rand01();
        double u2 = rand01();
        wi = sample_cosine_hemisphere(n, u1, u2);
    }
    // ---- Phong specular branch（只有在開啟時才有機會走到）----
#if ENABLE_PHONG_SPECULAR
    else {
        vec3 rdir = reflect_vec(-unit_vector(wo), n);
        if (rdir.length_squared() == 0.0) return 0;
        double u1 = rand01();
        double u2 = rand01();
        wi = sample_phong_lobe(rdir, mat->shininess, u1, u2);
        if (dot(n, wi) <= 0.0) wi = -wi;
    }
#endif

    double pdf_diff = cosine_hemisphere_pdf(n, wi);
    vec3 rdir = reflect_vec(-unit_vector(wo), n);
    double pdf_spec = 0.0;
#if ENABLE_PHONG_SPECULAR
    if (ks > 0.0) {
        pdf_spec = phong_lobe_pdf(mat->shininess, rdir, wi);
    }
#endif

    double pdf = pd * pdf_diff + ps * pdf_spec;
    if (pdf <= 0.0) return 0;

    *wi_out  = wi;
    *pdf_out = pdf;
    return 1;
}

#endif // TRACER_BSDF_SAMPLING_H
