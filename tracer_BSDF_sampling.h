#ifndef TRACER_BSDF_SAMPLING_H
#define TRACER_BSDF_SAMPLING_H

#include "Speed_up.h"
#include "tracer_light_source.h"  

//  LOCAL 座標轉換 

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

//  BSDF sampling（半球 / Phong lobe）  

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

//  BRDF eval（Lambert + Phong） 

inline vec3 eval_brdf(const Material* mat,
                      const vec3& n,
                      const vec3& wi,
                      const vec3& wo) {
    vec3 wi_n = unit_vector(wi);
    vec3 wo_n = unit_vector(wo);

    double cos_theta_i = dot(n, wi_n);
    if (cos_theta_i < 0.0) cos_theta_i = 0.0;
    if (cos_theta_i <= 0.0) return vec3(0.0, 0.0, 0.0);

    // Diffuse: ρ_d / π
    vec3 fd = (mat->kd / M_PI) * mat->color;

    // Specular (Phong)
    vec3 rdir = reflect_vec(-wo_n, n);
    double cos_alpha = dot(unit_vector(rdir), wi_n);
    if (cos_alpha < 0.0) cos_alpha = 0.0;
    vec3 fs(0.0,0.0,0.0);
    if (cos_alpha > 0.0 && mat->ks > 0.0) {
        double norm = (mat->shininess + 2.0) / (2.0 * M_PI);
        fs = mat->ks * norm * pow(cos_alpha, mat->shininess) * mat->color;
    }

    return fd + fs;
}

//  BRDF direction sampling 

inline int sample_brdf_direction(const Material* mat,
                                 const vec3& n,
                                 const vec3& wo,
                                 vec3* wi_out,
                                 double* pdf_out) {
    double kd = mat->kd;
    double ks = mat->ks;
    if (kd < 0.0) kd = 0.0;
    if (ks < 0.0) ks = 0.0;
    double sum = kd + ks;
    double pd = 1.0, ps = 0.0;
    if (sum > 0.0) {
        pd = kd / sum;
        ps = ks / sum;
    }

    double u = rand01();
    int use_diffuse = (u < pd);

    vec3 wi;
    if (use_diffuse) {
        double u1 = rand01();
        double u2 = rand01();
        wi = sample_cosine_hemisphere(n, u1, u2);
    } else {
        vec3 rdir = reflect_vec(-unit_vector(wo), n);
        if (rdir.length_squared() == 0.0) return 0;
        double u1 = rand01();
        double u2 = rand01();
        wi = sample_phong_lobe(rdir, mat->shininess, u1, u2);
        if (dot(n, wi) <= 0.0) wi = -wi;
    }

    double pdf_diff = cosine_hemisphere_pdf(n, wi);
    vec3 rdir = reflect_vec(-unit_vector(wo), n);
    double pdf_spec = phong_lobe_pdf(mat->shininess, rdir, wi);
    double pdf = pd * pdf_diff + ps * pdf_spec;
    if (pdf <= 0.0) return 0;

    *wi_out  = wi;
    *pdf_out = pdf;
    return 1;
}

#endif // TRACER_BSDF_SAMPLING_H
