// Last confirm date 2025/12/06 
#ifndef SPEED_UP_H
#define SPEED_UP_H

#include <float.h>
#include <math.h>
#include "vec3.h"
#include "tracer_camera.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef vec3 color;

//  Material / Geometry 

struct Material {
    vec3   color;
    double ka;
    double kd;
    double ks;
    double shininess;
    double reflectivity;
};

struct SphereLight {
    point3  center;
    double  radius;
    vec3    emission;
    Material mat;
};

struct Sphere {
    point3  center;
    double  radius;
    Material mat;
};

struct Triangle {
    point3 v0, v1, v2;
    vec3   normal;
    Material mat;
};

struct HitRecord {
    point3 p;
    vec3   normal;
    double t;
    int    front_face;
    const Material* mat;
    vec3   emission;
    int    is_light;
    int    light_index;
};

inline void set_face_normal(const Ray* r, const vec3& outward, HitRecord* rec) {
    rec->front_face = dot(r->dir, outward) < 0;
    rec->normal     = rec->front_face ? outward : -outward;
}

//  Scene + BVH 

#define MAX_LIGHTS 64
#define MAX_TRIS   2000500  // [修改] 加大到 10 萬，確保可以塞入超多細節
#define MAX_SPHERES  128

struct AABB {
    vec3 minimum;
    vec3 maximum;
};

inline double get_coord(const vec3& v, int axis) {
    if (axis == 0) return v.x();
    if (axis == 1) return v.y();
    return v.z();
}

inline AABB make_empty_box() {
    double inf = DBL_MAX;
    AABB b;
    b.minimum = vec3( inf,  inf,  inf);
    b.maximum = vec3(-inf, -inf, -inf);
    return b;
}

inline AABB expand_to_include(const AABB& box, const vec3& p) {
    AABB r;
    r.minimum = vec3(fmin(box.minimum.x(), p.x()),
                     fmin(box.minimum.y(), p.y()),
                     fmin(box.minimum.z(), p.z()));
    r.maximum = vec3(fmax(box.maximum.x(), p.x()),
                     fmax(box.maximum.y(), p.y()),
                     fmax(box.maximum.z(), p.z()));
    return r;
}

inline AABB surrounding_box(const AABB& a, const AABB& b) {
    AABB r;
    r.minimum = vec3(fmin(a.minimum.x(), b.minimum.x()),
                     fmin(a.minimum.y(), b.minimum.y()),
                     fmin(a.minimum.z(), b.minimum.z()));
    r.maximum = vec3(fmax(a.maximum.x(), b.maximum.x()),
                     fmax(a.maximum.y(), b.maximum.y()),
                     fmax(a.maximum.z(), b.maximum.z()));
    return r;
}

inline AABB triangle_bounds(const Triangle* tri) {
    AABB b = make_empty_box();
    b = expand_to_include(b, tri->v0);
    b = expand_to_include(b, tri->v1);
    b = expand_to_include(b, tri->v2);
    return b;
}

inline int aabb_hit(const AABB& box, const Ray* r, double t_min, double t_max) {
    for (int axis = 0; axis < 3; ++axis) {
        double invD = 1.0 / get_coord(r->dir, axis);
        double t0 = (get_coord(box.minimum, axis) - get_coord(r->orig, axis)) * invD;
        double t1 = (get_coord(box.maximum, axis) - get_coord(r->orig, axis)) * invD;
        if (invD < 0.0) {
            double tmp = t0; t0 = t1; t1 = tmp;
        }
        if (t0 > t_min) t_min = t0;
        if (t1 < t_max) t_max = t1;
        if (t_max <= t_min) return 0;
    }
    return 1;
}

struct BVHNode {
    AABB box;
    int  left;
    int  right;
    int  start;
    int  count;
    int  is_leaf;
};

struct Scene {
    SphereLight lights[MAX_LIGHTS];
    int         num_lights;

    Sphere      spheres[MAX_SPHERES];
    int         num_spheres;

    Triangle    triangles[MAX_TRIS];
    int         num_tris;

    // BVH over triangles
    int      tri_indices[MAX_TRIS];
    BVHNode  bvh_nodes[2 * MAX_TRIS];
    int      num_bvh_nodes;
};

//  Basic hit functions (per-primitive) 

inline int hit_sphere(const SphereLight* s, int light_index,
                      const Ray* r, double t_min, double t_max, HitRecord* rec) {
    vec3 oc = r->orig - s->center;
    double a = dot(r->dir, r->dir);
    double half_b = dot(oc, r->dir);
    double c = dot(oc, oc) - s->radius * s->radius;
    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return 0;
    double sqrtd = sqrt(discriminant);

    double root = (-half_b - sqrtd) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || root > t_max) return 0;
    }

    rec->t = root;
    rec->p = ray_at(r, rec->t);
    vec3 outward = (rec->p - s->center) / s->radius;
    set_face_normal(r, outward, rec);
    rec->mat       = &s->mat;
    rec->emission  = s->emission;
    rec->is_light  = (s->emission.x() > 0.0 || s->emission.y() > 0.0 || s->emission.z() > 0.0);
    rec->light_index = light_index;
    return 1;
}

// 一般非發光球的 intersection
inline int hit_sphere_geom(const Sphere* s,
                           const Ray* r, double t_min, double t_max,
                           HitRecord* rec) {
    vec3 oc = r->orig - s->center;
    double a = dot(r->dir, r->dir);
    double half_b = dot(oc, r->dir);
    double c = dot(oc, oc) - s->radius * s->radius;
    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return 0;
    double sqrtd = sqrt(discriminant);

    double root = (-half_b - sqrtd) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || root > t_max) return 0;
    }

    rec->t = root;
    rec->p = ray_at(r, rec->t);
    vec3 outward = (rec->p - s->center) / s->radius;
    set_face_normal(r, outward, rec);
    rec->mat       = &s->mat;
    rec->emission  = vec3(0.0, 0.0, 0.0);
    rec->is_light  = 0;
    rec->light_index = -1;
    return 1;
}

inline int hit_triangle(const Triangle* tri,
                        const Ray* r, double t_min, double t_max, HitRecord* rec) {
    const double EPS = 1e-6;
    vec3 v0v1 = tri->v1 - tri->v0;
    vec3 v0v2 = tri->v2 - tri->v0;
    vec3 pvec = cross(r->dir, v0v2);
    double det = dot(v0v1, pvec);

    if (fabs(det) < EPS) return 0;
    double invDet = 1.0 / det;

    vec3 tvec = r->orig - tri->v0;
    double u = dot(tvec, pvec) * invDet;
    if (u < 0.0 || u > 1.0) return 0;

    vec3 qvec = cross(tvec, v0v1);
    double v = dot(r->dir, qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) return 0;

    double t = dot(v0v2, qvec) * invDet;
    if (t < t_min || t > t_max) return 0;

    rec->t = t;
    rec->p = ray_at(r, rec->t);
    vec3 outward = tri->normal;
    set_face_normal(r, outward, rec);
    rec->mat       = &tri->mat;
    rec->emission  = vec3(0.0, 0.0, 0.0);
    rec->is_light  = 0;
    rec->light_index = -1;
    return 1;
}

//  BVH build & traversal 

static inline int bvh_build_node(Scene* scene, int start, int end) {
    int node_index = scene->num_bvh_nodes++;
    BVHNode* node = &scene->bvh_nodes[node_index];
    node->start = start;
    node->count = end - start;
    node->left  = -1;
    node->right = -1;
    node->is_leaf = 0;

    AABB bounds = make_empty_box();
    AABB centroid_bounds = make_empty_box();

    for (int i = start; i < end; ++i) {
        int tri_idx = scene->tri_indices[i];
        const Triangle* tri = &scene->triangles[tri_idx];
        AABB tb = triangle_bounds(tri);
        bounds = surrounding_box(bounds, tb);

        vec3 centroid = (tri->v0 + tri->v1 + tri->v2) / 3.0;
        centroid_bounds = expand_to_include(centroid_bounds, centroid);
    }
    node->box = bounds;

    int n = end - start;
    if (n <= 4) {
        node->is_leaf = 1;
        return node_index;
    }

    vec3 diag = centroid_bounds.maximum - centroid_bounds.minimum;
    int axis = 0;
    if (diag.y() > diag.x() && diag.y() >= diag.z()) axis = 1;
    else if (diag.z() > diag.x() && diag.z() > diag.y()) axis = 2;

    double mid = 0.5 * (get_coord(centroid_bounds.minimum, axis) +
                        get_coord(centroid_bounds.maximum, axis));

    int mid_index = start;
    for (int i = start; i < end; ++i) {
        int tri_idx = scene->tri_indices[i];
        const Triangle* tri = &scene->triangles[tri_idx];
        vec3 centroid = (tri->v0 + tri->v1 + tri->v2) / 3.0;
        if (get_coord(centroid, axis) < mid) {
            int tmp = scene->tri_indices[i];
            scene->tri_indices[i] = scene->tri_indices[mid_index];
            scene->tri_indices[mid_index] = tmp;
            mid_index++;
        }
    }

    if (mid_index == start || mid_index == end) {
        node->is_leaf = 1;
        return node_index;
    }

    node->left  = bvh_build_node(scene, start, mid_index);
    node->right = bvh_build_node(scene, mid_index, end);
    node->is_leaf = 0;
    return node_index;
}

inline void build_scene_bvh(Scene* scene) {
    scene->num_bvh_nodes = 0;
    if (scene->num_tris <= 0) return;
    for (int i = 0; i < scene->num_tris; ++i) {
        scene->tri_indices[i] = i;
    }
    bvh_build_node(scene, 0, scene->num_tris);
}

static inline int bvh_hit_node(const Scene* scene,
                               int node_index,
                               const Ray* r,
                               double t_min,
                               double t_max,
                               HitRecord* out_rec,
                               double* closest_t) {
    const BVHNode* node = &scene->bvh_nodes[node_index];
    if (!aabb_hit(node->box, r, t_min, t_max)) return 0;

    int hit_anything = 0;
    HitRecord temp;

    if (node->is_leaf) {
        for (int i = 0; i < node->count; ++i) {
            int tri_index = scene->tri_indices[node->start + i];
            const Triangle* tri = &scene->triangles[tri_index];
            if (hit_triangle(tri, r, t_min, *closest_t, &temp)) {
                hit_anything = 1;
                *closest_t = temp.t;
                *out_rec = temp;
            }
        }
    } else {
        if (node->left  >= 0)
            hit_anything |= bvh_hit_node(scene, node->left,  r, t_min, *closest_t, out_rec, closest_t);
        if (node->right >= 0)
            hit_anything |= bvh_hit_node(scene, node->right, r, t_min, *closest_t, out_rec, closest_t);
    }
    return hit_anything;
}

//  Scene hit using BVH + lights 

inline int hit_scene(const Scene* scene,
                     const Ray* r, double t_min, double t_max, HitRecord* rec) {
    HitRecord temp;
    int hit_anything = 0;
    double closest = t_max;

    // Triangles via BVH
    if (scene->num_bvh_nodes > 0) {
        if (bvh_hit_node(scene, 0, r, t_min, closest, &temp, &closest)) {
            hit_anything = 1;
            *rec = temp;
        }
    } else {
        for (int i = 0; i < scene->num_tris; ++i) {
            if (hit_triangle(&scene->triangles[i], r, t_min, closest, &temp)) {
                hit_anything = 1;
                closest = temp.t;
                *rec = temp;
            }
        }
    }

    // 一般非發光球
    for (int i = 0; i < scene->num_spheres; ++i) {
        if (hit_sphere_geom(&scene->spheres[i], r, t_min, closest, &temp)) {
            hit_anything = 1;
            closest = temp.t;
            *rec = temp;
        }
    }

    // Lights (emissive spheres)
    for (int i = 0; i < scene->num_lights; ++i) {
        if (hit_sphere(&scene->lights[i], i, r, t_min, closest, &temp)) {
            hit_anything = 1;
            closest = temp.t;
            *rec = temp;
        }
    }

    return hit_anything;
}

#endif // SPEED_UP_H
