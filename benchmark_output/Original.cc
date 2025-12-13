// Last confirm date 2025/12/06 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "vec3.h"
#include "tracer_camera.h"
#include "Speed_up.h"
#include "tracer_MIS.h"

#include <mach-o/dyld.h>
#include <limits.h>
#include <unistd.h>
#include <libgen.h>


#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef OUT_DIR
#define OUT_DIR ""
#endif

//  Scene loading: E / V / F / R / M / SL / T 

int load_scene_from_file(const char* filename,
                         Scene* scene,
                         Camera* cam,
                         int* out_width,
                         int* out_height) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open input file: %s\n", filename);
        return 0;
    }

    scene->num_lights    = 0;
    scene->num_spheres   = 0;
    scene->num_tris      = 0;
    scene->num_bvh_nodes = 0;

    Material currentMat;
    currentMat.color        = vec3(1.0,1.0,1.0);
    //currentMat.ka           = 0.0;
    currentMat.kd           = 1.0;
    currentMat.ks           = 0.0;
    currentMat.shininess    = 16.0;
    currentMat.reflectivity = 0.0;

    point3 eye;
    vec3   view_dir;
    vec3   up_vec;
    double fov_deg = 0.0;
    int    file_width  = -1;
    int    file_height = -1;

    int gotE = 0;
    int gotV = 0;
    int gotF = 0;
    int gotR = 0;

    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        char* hash = strchr(line, '#');
        if (hash) *hash = '\0';

        char* p = line;
        while (*p==' ' || *p=='\t' || *p=='\r' || *p=='\n') ++p;
        if (*p == '\0') continue;

        char tag[8];
        if (sscanf(p, "%7s", tag) != 1) continue;

        if (strcmp(tag,"E") == 0) {
            double ex,ey,ez;
            if (sscanf(p, "%*s %lf %lf %lf", &ex,&ey,&ez) == 3) {
                eye = point3(ex,ey,ez);
                gotE = 1;
            } else {
                fprintf(stderr, "Error: bad 'E' line in %s\n", filename);
            }
        } else if (strcmp(tag,"V") == 0) {
            double vx,vy,vz, ux,uy,uz;
            if (sscanf(p, "%*s %lf %lf %lf %lf %lf %lf",
                       &vx,&vy,&vz,&ux,&uy,&uz) == 6) {
                view_dir = vec3(vx,vy,vz);
                up_vec   = vec3(ux,uy,uz);
                gotV = 1;
            } else {
                fprintf(stderr, "Error: bad 'V' line in %s\n", filename);
            }
        } else if (strcmp(tag,"F") == 0) {
            double fdeg;
            if (sscanf(p, "%*s %lf", &fdeg) == 1) {
                fov_deg = fdeg;
                gotF = 1;
            } else {
                fprintf(stderr, "Error: bad 'F' line in %s\n", filename);
            }
        } else if (strcmp(tag,"R") == 0) {
            int w,h;
            if (sscanf(p, "%*s %d %d", &w,&h) == 2 && w > 0 && h > 0) {
                file_width  = w;
                file_height = h;
                gotR = 1;
            } else {
                fprintf(stderr, "Error: bad 'R' line in %s\n", filename);
            }
        } else if (strcmp(tag,"M") == 0) {
            double cr, cg, cb, kd, ks, shininess, refl;
            if (sscanf(p, "%*s %lf %lf %lf %lf %lf %lf %lf",
                       &cr,&cg,&cb,&kd,&ks,&shininess,&refl) == 7) {
                currentMat.color        = vec3(cr,cg,cb);
                //currentMat.ka           = ka;
                currentMat.kd           = kd;
                currentMat.ks           = ks;
                currentMat.shininess    = shininess;
                currentMat.reflectivity = refl;
            } else {
                fprintf(stderr, "Warning: bad 'M' line in %s, ignoring.\n", filename);
            }
        } else if (strcmp(tag,"SL") == 0) {
            double x,y,z,radius, er,eg,eb;
            if (sscanf(p, "%*s %lf %lf %lf %lf %lf %lf %lf",
                       &x,&y,&z,&radius,&er,&eg,&eb) == 7) {
                if (scene->num_lights < MAX_LIGHTS) {       //max light = 64
                    SphereLight* L = &scene->lights[scene->num_lights++];
                    L->center   = point3(x,y,z);
                    L->radius   = radius;
                    L->emission = vec3(er,eg,eb);
                    L->mat      = currentMat;
                } else {
                    fprintf(stderr, "Warning: too many lights, ignoring extra.\n");
                }
            } else {
                fprintf(stderr, "Warning: bad 'SL' line in %s, ignoring.\n", filename);
            }
        } else if (strcmp(tag,"S") == 0) {
            double ox,oy,oz,radius;
            if (sscanf(p, "%*s %lf %lf %lf %lf",
                       &ox,&oy,&oz,&radius) == 4) {
                if (scene->num_spheres < MAX_SPHERES) {
                    Sphere* s = &scene->spheres[scene->num_spheres++];
                    s->center = point3(ox,oy,oz);
                    s->radius = radius;
                    s->mat    = currentMat;
                } else {
                    fprintf(stderr, "Warning: too many spheres, ignoring extra.\n");
                }
            } else {
                fprintf(stderr, "Warning: bad 'S' line in %s, ignoring.\n", filename);
            }

        } else if (strcmp(tag,"T") == 0) {      //max triangles = 64
            double x0,y0,z0,x1,y1,z1,x2,y2,z2;
            if (sscanf(p, "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       &x0,&y0,&z0,&x1,&y1,&z1,&x2,&y2,&z2) == 9) {
                if (scene->num_tris < MAX_TRIS) {
                    Triangle* tri = &scene->triangles[scene->num_tris++];
                    tri->v0 = point3(x0,y0,z0);
                    tri->v1 = point3(x1,y1,z1);
                    tri->v2 = point3(x2,y2,z2);
                    tri->normal = unit_vector(cross(tri->v1 - tri->v0,
                                                    tri->v2 - tri->v0));
                    tri->mat    = currentMat;
                } else {
                    fprintf(stderr, "Warning: too many triangles, ignoring extra.\n");
                }
            } else {
                fprintf(stderr, "Warning: bad 'T' line in %s, ignoring.\n", filename);
            }
        }
    }

    fclose(fp);

    // 檢查 camera 必要項目
    int missing = 0;
    if (!gotE) {
        fprintf(stderr, "Error: missing 'E' (eye position) in scene file %s\n", filename);
        missing = 1;
    }
    if (!gotV) {
        fprintf(stderr, "Error: missing 'V' (view direction / up vector) in scene file %s\n", filename);
        missing = 1;
    }
    if (!gotF) {
        fprintf(stderr, "Error: missing 'F' (field-of-view) in scene file %s\n", filename);
        missing = 1;
    }
    if (!gotR) {
        fprintf(stderr, "Error: missing 'R' (image resolution) in scene file %s\n", filename);
        missing = 1;
    }
    if (missing) {
        return 0;
    }

    // 建 BVH
    build_scene_bvh(scene);

    // 使用 input 檔的參數建立 Camera
    setup_camera(cam, eye, view_dir, up_vec, fov_deg, file_width, file_height);
    *out_width  = file_width;
    *out_height = file_height;

    return 1;
}

//  Rendering 

void render_image(const Scene* scene,
                  const Camera* cam,
                  const char* filename,
                  RenderMode mode,
                  int image_width,
                  int image_height,
                  int samples_per_pixel,
                  int max_depth) {
    FILE* out = fopen(filename, "w");
    if (!out) {
        fprintf(stderr, "Failed to open output file: %s\n", filename);
        return;
    }

    const int MAX_COLOR = 256; // 依你原本指定 256

    fprintf(out, "P3\n%d %d\n%d\n", image_width, image_height, MAX_COLOR);

    for (int j = image_height - 1; j >= 0; --j) {
        fprintf(stderr, "\rRendering %s line %d/%d",
                filename, image_height - 1 - j, image_height);
        fflush(stderr);
        for (int i = 0; i < image_width; ++i) {
            vec3 col(0.0,0.0,0.0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                double u = (i + rand01()) / (double)image_width;
                double v = (j + rand01()) / (double)image_height;
                Ray r = get_camera_ray(cam, u, v);
                col += ray_color(scene, &r, mode, max_depth);
            }
            col /= (double)samples_per_pixel;
            col = tonemap(col);
            col = gamma_correct(col, 2.2);

            double rr = col.x();
            double gg = col.y();
            double bb = col.z();
            if (rr < 0.0) rr = 0.0; else if (rr > 1.0) rr = 1.0;
            if (gg < 0.0) gg = 0.0; else if (gg > 1.0) gg = 1.0;
            if (bb < 0.0) bb = 0.0; else if (bb > 1.0) bb = 1.0;

            int ir = (int)(rr * MAX_COLOR + 0.5);
            int ig = (int)(gg * MAX_COLOR + 0.5);
            int ib = (int)(bb * MAX_COLOR + 0.5);

            fprintf(out, "%d %d %d\n", ir, ig, ib);
        }
    }
    fprintf(stderr, "\nDone %s\n", filename);
    fclose(out);
}

// G: light sampling; Red: BSDF
void render_mis_weight_image(const Scene* scene,
                             const Camera* cam,
                             const char* filename,
                             int image_width,
                             int image_height,
                             int samples_per_pixel) {
    FILE* out = fopen(filename, "w");
    if (!out) {
        fprintf(stderr, "Failed to open output file: %s\n", filename);
        return;
    }

    const int MAX_COLOR = 255;

    fprintf(out, "P3\n%d %d\n%d\n", image_width, image_height, MAX_COLOR);

    for (int j = image_height - 1; j >= 0; --j) {
        fprintf(stderr, "\rRendering %s (MIS weight) line %d/%d",
                filename, image_height - 1 - j, image_height);
        fflush(stderr);

        for (int i = 0; i < image_width; ++i) {
            vec3 col(0.0,0.0,0.0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                double u = (i + rand01()) / (double)image_width;
                double v = (j + rand01()) / (double)image_height;
                Ray r = get_camera_ray(cam, u, v);
                col += ray_mis_weight_color(scene, &r, 1);
            }
            col /= (double)samples_per_pixel;

            double rr = col.x();
            double gg = col.y();
            double bb = col.z();
            if (rr < 0.0) rr = 0.0; else if (rr > 1.0) rr = 1.0;
            if (gg < 0.0) gg = 0.0; else if (gg > 1.0) gg = 1.0;
            if (bb < 0.0) bb = 0.0; else if (bb > 1.0) bb = 1.0;

            int ir = (int)(rr * MAX_COLOR + 0.5);
            int ig = (int)(gg * MAX_COLOR + 0.5);
            int ib = (int)(bb * MAX_COLOR + 0.5);

            fprintf(out, "%d %d %d\n", ir, ig, ib);
        }
    }

    fprintf(stderr, "\nDone %s\n", filename);
    fclose(out);
}


//  main 

int main(int argc, char** argv) {
    // 保證 Finder 雙擊時，相對路徑會以執行檔所在目錄為基準
    {
        char exe_path[PATH_MAX];
        uint32_t size = sizeof(exe_path);

        if (_NSGetExecutablePath(exe_path, &size) == 0) {
            // chdir(dirname(exe_path));
        }
    }


    const char* input_file = "input.txt";
    if (argc > 1) {
        input_file = argv[1];
    }

    static Scene scene;

    Camera cam;
    int image_width  = 0;
    int image_height = 0;

    if (!load_scene_from_file(input_file, &scene, &cam,
                              &image_width, &image_height)) {
        fprintf(stderr, "Scene loading failed.\n");
        return 1;
    }

    srand(42);

    const int SAMPLES_PER_PIXEL  = 30;
    const int MAX_DEPTH          = 5;

    /* render_image(&scene, &cam, OUT_DIR "AdvCG_light.ppm",
                 MODE_LIGHT,
                 image_width, image_height,
                 SAMPLES_PER_PIXEL, MAX_DEPTH); */

    /* render_image(&scene, &cam, OUT_DIR "AdvCG_bsdf.ppm",
                 MODE_BRDF,
                 image_width, image_height,
                 SAMPLES_PER_PIXEL, MAX_DEPTH); */

    render_image(&scene, &cam, OUT_DIR "AdvCG_mis.ppm",
                 MODE_MIS,
                 image_width, image_height,
                 SAMPLES_PER_PIXEL, MAX_DEPTH);

    /* render_mis_weight_image(&scene, &cam, OUT_DIR "AdvCG_mis_weight.ppm",
                            image_width, image_height,
                            SAMPLES_PER_PIXEL); */    

    return 0;
}


