// Last confirm date 2025/12/06 
#ifndef TRACER_MIS_H
#define TRACER_MIS_H

#include "tracer_BSDF_sampling.h" // 內含 Speed_up.h + tracer_light_source.h

//  MIS 資料結構 

struct MISSample {
    vec3 wi;        // 入射方向
    vec3 Li;        // 來自光源的 radiance
    vec3 f;         // BRDF 評估值
    double cos_theta;
    double pdf_light;
    double pdf_brdf;
};

inline double balance_weight(double n_a, double pdf_a,
                             double n_b, double pdf_b) {
    double a = n_a * pdf_a;
    double b = n_b * pdf_b;
    double denom = a + b;
    if (denom <= 0.0) return 0.0;
    return a / denom;
}

// evaluate_mis：
//   - 回傳最終 radiance L
//   - 若 ratio_color != 0，順便算出：
//       R = BRDF-sampling 貢獻比例
//       G = light-sampling 貢獻比例
inline vec3 evaluate_mis(const MISSample* light_samples, int n_light_samples,
                         const MISSample* brdf_samples,  int n_brdf_samples,
                         double n_light, double n_brdf,
                         vec3* ratio_color = 0) {
    vec3 L(0.0,0.0,0.0);

    vec3 L_from_light(0.0,0.0,0.0);
    vec3 L_from_brdf (0.0,0.0,0.0);

    // Light-sampling branch
    for (int i = 0; i < n_light_samples; ++i) {
        const MISSample* s = &light_samples[i];
        if (s->pdf_light <= 0.0 || s->cos_theta <= 0.0) continue;

        double w = balance_weight(n_light, s->pdf_light,
                                  n_brdf,  s->pdf_brdf);

        vec3 contrib = s->f * s->Li * (s->cos_theta * w / s->pdf_light);
        L           += contrib;
        L_from_light += contrib;
    }

    // BRDF-sampling branch
    for (int i = 0; i < n_brdf_samples; ++i) {
        const MISSample* s = &brdf_samples[i];
        if (s->pdf_brdf <= 0.0 || s->cos_theta <= 0.0) continue;

        double w = balance_weight(n_brdf,  s->pdf_brdf,
                                  n_light, s->pdf_light);

        vec3 contrib = s->f * s->Li * (s->cos_theta * w / s->pdf_brdf);
        L           += contrib;
        L_from_brdf += contrib;
    }

    if (ratio_color) {
        // 用簡單的 luminance 近似衡量能量大小
        auto luminance = [](const vec3& c) {
            return 0.2126 * c.x() + 0.7152 * c.y() + 0.0722 * c.z();
        };
        double lum_light = luminance(L_from_light);
        double lum_brdf  = luminance(L_from_brdf);
        double lum_sum   = lum_light + lum_brdf;

        if (lum_sum > 0.0) {
            double r = lum_brdf  / lum_sum;  // BRDF 比例
            double g = lum_light / lum_sum;  // Light 比例
            *ratio_color = vec3(r, g, 0.0);
        } else {
            *ratio_color = vec3(0.0,0.0,0.0);
        }
    }

    return L;
}

//  Render 模式 

enum RenderMode {
    MODE_LIGHT = 0, // 只 light sampling
    MODE_BRDF  = 1, // 只 BRDF sampling
    MODE_MIS   = 2  // 兩者 MIS
};

#define MAX_LIGHT_SAMPLES 64
#define MAX_BRDF_SAMPLES  16

//  共用：建立 MIS sample 集合 

inline void gather_mis_samples(const Scene* scene,
                               const HitRecord* rec,
                               const vec3& wo,
                               RenderMode mode,
                               MISSample* light_samples,
                               int* n_light_samples,
                               MISSample* brdf_samples,
                               int* n_brdf_samples,
                               double* n_light,
                               double* n_brdf) {
    *n_light_samples = 0;
    *n_brdf_samples  = 0;

    int Nl = (mode == MODE_BRDF) ? 0 : 1;  // 有沒有做 light sampling
    int Nb = (mode == MODE_LIGHT) ? 0 : 1; // 有沒有做 BRDF sampling

    *n_light = (double)Nl;
    *n_brdf  = (double)Nb;

    // ---------- Light sampling ----------
    if (Nl > 0) {
        for (int i = 0; i < scene->num_lights && *n_light_samples < MAX_LIGHT_SAMPLES; ++i) {
            const SphereLight* L = &scene->lights[i];
            vec3 wi;
            double dist;
            double pdf_light;
            vec3 Li;
            if (!sample_sphere_light_from_point(rec->p, rec->normal, L,
                                                &wi, &dist, &pdf_light, &Li)) {
                continue;
            }
            if (pdf_light <= 0.0) continue;
            if (scene_occluded(scene, rec->p, wi, dist)) continue;

            double cos_theta = dot(rec->normal, wi);
            if (cos_theta < 0.0) cos_theta = 0.0;
            if (cos_theta <= 0.0) continue;

            vec3 f = eval_brdf(rec->mat, rec->normal, wi, wo);

            double kd, ks;
            get_effective_kd_ks(rec->mat, &kd, &ks);
            double sum = kd + ks;
            double pd = 1.0, ps = 0.0;
            if (sum > 0.0) {
                pd = kd / sum;
                ps = ks / sum;
            }
            double pdf_diff = cosine_hemisphere_pdf(rec->normal, wi);
            vec3 rdir = reflect_vec(-unit_vector(wo), rec->normal);
            double pdf_spec = phong_lobe_pdf(rec->mat->shininess, rdir, wi);
            double pdf_brdf = pd * pdf_diff + ps * pdf_spec;

            MISSample s;
            s.wi = wi;
            s.Li = Li;
            s.f  = f;
            s.cos_theta = cos_theta;
            s.pdf_light = pdf_light;
            s.pdf_brdf  = pdf_brdf;

            light_samples[(*n_light_samples)++] = s;
        }
    }

    // ---------- BRDF sampling ----------
    if (Nb > 0 && *n_brdf_samples < MAX_BRDF_SAMPLES) {
        vec3 wi;
        double pdf_brdf;
        if (sample_brdf_direction(rec->mat, rec->normal, wo, &wi, &pdf_brdf)) {
            if (pdf_brdf > 0.0 && dot(rec->normal, wi) > 0.0) {
                vec3 Li;
                double pdf_light;
                if (trace_to_light(scene, rec->p, wi, &Li, &pdf_light)) {
                    double cos_theta = dot(rec->normal, wi);
                    if (cos_theta < 0.0) cos_theta = 0.0;
                    vec3 f = eval_brdf(rec->mat, rec->normal, wi, wo);

                    MISSample s;
                    s.wi = wi;
                    s.Li = Li;
                    s.f  = f;
                    s.cos_theta = cos_theta;
                    s.pdf_light = pdf_light;
                    s.pdf_brdf  = pdf_brdf;

                    brdf_samples[(*n_brdf_samples)++] = s;
                }
            }
        }
    }
}

//  Direct lighting（正常 shading 用） 

inline vec3 direct_lighting(const Scene* scene,
                            const HitRecord* rec,
                            const vec3& wo,
                            RenderMode mode) {
    MISSample light_samples[MAX_LIGHT_SAMPLES];
    MISSample brdf_samples[MAX_BRDF_SAMPLES];
    int n_light_samples = 0;
    int n_brdf_samples  = 0;
    double n_light = 0.0;
    double n_brdf  = 0.0;

    gather_mis_samples(scene, rec, wo, mode,
                       light_samples, &n_light_samples,
                       brdf_samples, &n_brdf_samples,
                       &n_light, &n_brdf);

    return evaluate_mis(light_samples, n_light_samples,
                        brdf_samples, n_brdf_samples,
                        n_light, n_brdf, 0);
}

//  Direct lighting（MIS 權重圖用） 
// 回傳顏色：R=BRDF-sampling 能量比例、G=light-sampling 能量比例
inline vec3 direct_lighting_mis_weight(const Scene* scene,
                                       const HitRecord* rec,
                                       const vec3& wo) {
    MISSample light_samples[MAX_LIGHT_SAMPLES];
    MISSample brdf_samples[MAX_BRDF_SAMPLES];
    int n_light_samples = 0;
    int n_brdf_samples  = 0;
    double n_light = 0.0;
    double n_brdf  = 0.0;

    // 權重圖固定看 MIS 模式：兩種 sampling 都啟用
    gather_mis_samples(scene, rec, wo, MODE_MIS,
                       light_samples, &n_light_samples,
                       brdf_samples, &n_brdf_samples,
                       &n_light, &n_brdf);

    vec3 ratio_color(0.0,0.0,0.0);
    (void) evaluate_mis(light_samples, n_light_samples,
                        brdf_samples, n_brdf_samples,
                        n_light, n_brdf,
                        &ratio_color);
    return ratio_color;
}

//  Tone mapping / Gamma / background 

inline vec3 tonemap(const vec3& c) {
    vec3 t;
    t.e[0] = c.e[0] / (1.0 + c.e[0]);
    t.e[1] = c.e[1] / (1.0 + c.e[1]);
    t.e[2] = c.e[2] / (1.0 + c.e[2]);
    return t;
}

inline vec3 gamma_correct(const vec3& c, double gamma) {
    double inv = 1.0 / gamma;
    return vec3(pow(fmax(0.0, c.e[0]), inv),
                pow(fmax(0.0, c.e[1]), inv),
                pow(fmax(0.0, c.e[2]), inv));
}

inline vec3 ray_background(const Ray* r) {
    (void)r;
    return vec3(0.0, 0.0, 0.0);
}

//  Path tracing: ray_color（正常渲染） 

inline vec3 ray_color(const Scene* scene,
                      const Ray* r,
                      RenderMode mode,
                      int depth) {
    const double EPS = 1e-4;
    if (depth <= 0) {
        return vec3(0.0,0.0,0.0);
    }

    HitRecord rec;
    if (!hit_scene(scene, r, EPS, DBL_MAX, &rec)) {
        return ray_background(r);
    }

    if (rec.is_light) {
        return rec.emission;
    }

    vec3 wo = -unit_vector(r->dir);

    vec3 ambient = rec.mat->ka * rec.mat->color;
    vec3 direct  = direct_lighting(scene, &rec, wo, mode);

    vec3 refl_color(0.0,0.0,0.0);
    if (rec.mat->reflectivity > 0.0) {
        vec3 refl_dir = reflect_vec(-unit_vector(r->dir), rec.normal);
        refl_dir = unit_vector(refl_dir);
        Ray refl_ray;
        refl_ray.orig = rec.p + EPS * refl_dir;
        refl_ray.dir  = refl_dir;
        refl_color = rec.mat->reflectivity *
                     ray_color(scene, &refl_ray, mode, depth - 1);
    }

    return ambient + direct + refl_color;
}

//  ray_mis_weight_color：只算 MIS 比例顏色 

inline vec3 ray_mis_weight_color(const Scene* scene,
                                 const Ray* r,
                                 int depth) {
    (void)depth; // 權重圖只看第一個交點的 direct lighting

    const double EPS = 1e-4;
    HitRecord rec;
    if (!hit_scene(scene, r, EPS, DBL_MAX, &rec)) {
        return vec3(0.0,0.0,0.0);
    }

    // 光源本身不畫權重（給 0）
    if (rec.is_light) {
        return vec3(0.0,0.0,0.0);
    }

    vec3 wo = -unit_vector(r->dir);
    return direct_lighting_mis_weight(scene, &rec, wo);
}

#endif // TRACER_MIS_H
