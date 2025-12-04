#ifndef TRACER_CAMERA_H
#define TRACER_CAMERA_H

#include <math.h>
#include "vec3.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Ray
struct Ray {
    point3 orig;
    vec3   dir;
};

inline point3 ray_at(const Ray* r, double t) {
    return r->orig + t * r->dir;
}

// Camera
struct Camera {
    point3 origin;
    vec3   forward;
    vec3   right;
    vec3   up;
    double viewport_width;
    double viewport_height;
};

// 依照 eye / view_dir / up_vec / fov / 解析度 建立 Camera
inline void setup_camera(Camera* cam,
                         const point3& eye,
                         const vec3& view_dir,
                         const vec3& up_vec,
                         double fov_deg,
                         int image_width,
                         int image_height) {
    double aspect_ratio = (double)image_width / (double)image_height;
    double fov_rad = fov_deg * M_PI / 180.0;
    double viewport_height = 2.0 * tan(0.5 * fov_rad);
    double viewport_width  = viewport_height * aspect_ratio;

    vec3 forward = unit_vector(view_dir);
    vec3 right   = unit_vector(cross(forward, up_vec));
    vec3 cam_up  = cross(right, forward);

    cam->origin  = eye;
    cam->forward = forward;
    cam->right   = right;
    cam->up      = cam_up;
    cam->viewport_width  = viewport_width;
    cam->viewport_height = viewport_height;
}

inline Ray get_camera_ray(const Camera* cam, double u, double v) {
    Ray r;
    vec3 dir = cam->forward
             + (2.0 * u - 1.0) * (cam->viewport_width  * 0.5) * cam->right
             + (2.0 * v - 1.0) * (cam->viewport_height * 0.5) * cam->up;
    r.orig = cam->origin;
    r.dir  = unit_vector(dir);
    return r;
}

#endif // TRACER_CAMERA_H
