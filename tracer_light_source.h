// Last confirm date 2025/12/06 
#ifndef TRACER_LIGHT_SOURCE_H
#define TRACER_LIGHT_SOURCE_H

#include <stdlib.h>
#include "Speed_up.h"

// Random
inline double rand01() {
    return (double)rand() / (double)RAND_MAX;
}

// Sphere light sampling
inline int sample_sphere_light_from_point(const point3& p, const vec3& n,
                                          const SphereLight* L,
                                          vec3* wi_out,
                                          double* dist_out,
                                          double* pdf_omega_out,
                                          vec3* Li_out) {
    (void)n;

    double u1 = rand01();
    double u2 = rand01();
    double z  = 1.0 - 2.0 * u1;
    double r  = sqrt(fmax(0.0, 1.0 - z*z));
    double phi = 2.0 * M_PI * u2;
    double x = r * cos(phi);
    double y = r * sin(phi);
    vec3 dir_surface(x,y,z);

    point3 x_L = L->center + L->radius * dir_surface;
    vec3 to_light = x_L - p;
    double dist = to_light.length();
    if (dist <= 0.0) return 0;
    vec3 wi = to_light / dist;

    vec3 n_L = unit_vector(x_L - L->center);
    double cos_theta_L = dot(n_L, -wi);
    if (cos_theta_L <= 0.0) return 0;

    double area = 4.0 * M_PI * L->radius * L->radius;
    double pdf_area = 1.0 / area;
    double pdf_omega = pdf_area * (dist * dist) / cos_theta_L;

    *wi_out = wi;
    *dist_out = dist;
    *pdf_omega_out = pdf_omega;
    *Li_out = L->emission;
    return 1;
}

// Occlusion & trace to light
inline int scene_occluded(const Scene* scene,
                          const point3& p,
                          const vec3& wi,
                          double max_dist) {
    const double EPS = 1e-4;
    Ray shadow_ray;
    shadow_ray.orig = p + EPS * wi;
    shadow_ray.dir  = wi;
    HitRecord rec;
    if (hit_scene(scene, &shadow_ray, EPS, max_dist - EPS, &rec)) {
        return 1;
    }
    return 0;
}

inline int trace_to_light(const Scene* scene,
                          const point3& p,
                          const vec3& wi,
                          vec3* Li_out,
                          double* pdf_light_out) {
    const double EPS = 1e-4;
    Ray r;
    r.orig = p + EPS * wi;
    r.dir  = wi;
    HitRecord rec;
    if (!hit_scene(scene, &r, EPS, DBL_MAX, &rec)) return 0;
    if (!rec.is_light || rec.light_index < 0) return 0;

    const SphereLight* L = &scene->lights[rec.light_index];

    point3 x_L = rec.p;
    vec3 to_light = x_L - p;
    double dist2 = to_light.length_squared();
    if (dist2 <= 0.0) return 0;
    vec3 dir = unit_vector(to_light);
    vec3 n_L = unit_vector(x_L - L->center);
    double cos_theta_L = dot(n_L, -dir);
    if (cos_theta_L <= 0.0) return 0;

    double area = 4.0 * M_PI * L->radius * L->radius;
    double pdf_omega = dist2 / (area * cos_theta_L);

    *Li_out = L->emission;
    *pdf_light_out = pdf_omega;
    return 1;
}

#endif // TRACER_LIGHT_SOURCE_H
