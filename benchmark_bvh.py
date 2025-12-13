import os
import time
import subprocess
import re
import math
import platform
import sys
import shutil

# ================= 參數設定 =================
OUTPUT_DIR = "benchmark_output"
COMPILER = "g++"
SOURCE_FILE = "AdvCG_Final_MIS.cc"
HEADER_FILE = "Speed_up.h"
SOURCE_INPUT = "input.txt" # 將自動搜尋
TARGET_INPUT = os.path.join(OUTPUT_DIR, "input_optimized.txt")

# [設定]
# 1. 解析度: 強制設為 1024x1024 以放大運算量，凸顯加速差異
TARGET_WIDTH = 2048         
TARGET_HEIGHT = 2048        

print(f"🔥 [Benchmark V17] Visual Identity (Spheres) + Shadow Cache + AnyHit + FastMath")
print(f"   Config: Resolution={TARGET_WIDTH}x{TARGET_HEIGHT}")

if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
exe_ext = ".exe" if platform.system() == "Windows" else ""
EXE_FAST = os.path.join(OUTPUT_DIR, f"render_original{exe_ext}")
EXE_OPT  = os.path.join(OUTPUT_DIR, f"render_optimized_v17{exe_ext}")

# ================= 1. 場景準備 (僅修改解析度) =================
def find_input_file():
    global SOURCE_INPUT
    candidates = ["input.txt", "Graphics-hw4_Final/input.txt"]
    for f in candidates:
        if os.path.exists(f):
            SOURCE_INPUT = f
            print(f"📄 Found input file at: {SOURCE_INPUT}")
            return True
    files = [f for f in os.listdir('.') if 'input' in f and f.endswith('.txt')]
    if files:
        SOURCE_INPUT = files[0]
        print(f"📄 Found input file at: {SOURCE_INPUT}")
        return True
    return False

def prepare_scene():
    if not find_input_file(): return False
    print(f"⚙️ [1/5] Preparing Scene (Resolution 1024x1024)...")
    
    with open(SOURCE_INPUT, 'r', encoding='utf-8', errors='ignore') as f_in, \
         open(TARGET_INPUT, 'w', encoding='utf-8') as f_out:
        
        for line in f_in:
            line_clean = line.strip()
            parts = line_clean.split()
            if not parts: 
                f_out.write(line)
                continue
            
            # 強制修改解析度
            if parts[0] == 'R':
                f_out.write(f"R {TARGET_WIDTH} {TARGET_HEIGHT}\n")
            else:
                f_out.write(line)
    
    print(f"   Done. Scene ready: {TARGET_INPUT}")
    return True

# ================= 2. 生成 Speed_up_opt.h (含 Shadow Cache) =================

CODE_SPEEDUP_OPT = r"""
#ifndef SPEED_UP_OPT_H
#define SPEED_UP_OPT_H

#include <algorithm>
#include <vector>
#include <float.h>
#include <math.h>
#include "vec3.h"
#include "tracer_camera.h"

#define MAX_LIGHTS 64
#define MAX_TRIS 2000500 
#define MAX_SPHERES 128

// Fast RNG
static unsigned int _rng_seed = 123456789;
inline double fast_rand01() {
    unsigned int x = _rng_seed; 
    x ^= x << 13; x ^= x >> 17; x ^= x << 5; 
    _rng_seed = x;
    return x * (1.0 / 4294967295.0);
}

struct Material { vec3 color; double ka; double kd; double ks; double shininess; double reflectivity; };
struct SphereLight { point3 center; double radius; vec3 emission; Material mat; };
struct Sphere { point3 center; double radius; Material mat; };
struct Triangle { point3 v0, v1, v2; vec3 normal; Material mat; };
struct HitRecord { point3 p; vec3 normal; double t; int front_face; const Material* mat; vec3 emission; int is_light; int light_index; };

inline void set_face_normal(const Ray* r, const vec3& outward, HitRecord* rec) {
    rec->front_face = dot(r->dir, outward) < 0;
    rec->normal = rec->front_face ? outward : -outward;
}

struct AABB { vec3 minimum; vec3 maximum; };
inline double get_coord(const vec3& v, int axis) { return (axis==0) ? v.x() : ((axis==1) ? v.y() : v.z()); }

inline AABB make_empty_box() {
    double inf = DBL_MAX;
    AABB b; b.minimum = vec3(inf,inf,inf); b.maximum = vec3(-inf,-inf,-inf);
    return b;
}
inline AABB expand_to_include(const AABB& box, const vec3& p) {
    AABB r;
    r.minimum = vec3(fmin(box.minimum.x(), p.x()), fmin(box.minimum.y(), p.y()), fmin(box.minimum.z(), p.z()));
    r.maximum = vec3(fmax(box.maximum.x(), p.x()), fmax(box.maximum.y(), p.y()), fmax(box.maximum.z(), p.z()));
    return r;
}
inline AABB surrounding_box(const AABB& a, const AABB& b) {
    AABB r;
    r.minimum = vec3(fmin(a.minimum.x(), b.minimum.x()), fmin(a.minimum.y(), b.minimum.y()), fmin(a.minimum.z(), b.minimum.z()));
    r.maximum = vec3(fmax(a.maximum.x(), b.maximum.x()), fmax(a.maximum.y(), b.maximum.y()), fmax(a.maximum.z(), b.maximum.z()));
    return r;
}
inline AABB triangle_bounds(const Triangle* tri) {
    AABB b = make_empty_box(); b = expand_to_include(b, tri->v0); b = expand_to_include(b, tri->v1); b = expand_to_include(b, tri->v2); return b;
}

inline int aabb_hit(const AABB& box, const Ray* r, double t_min, double t_max) {
    for (int axis = 0; axis < 3; ++axis) {
        double invD = 1.0 / get_coord(r->dir, axis);
        double t0 = (get_coord(box.minimum, axis) - get_coord(r->orig, axis)) * invD;
        double t1 = (get_coord(box.maximum, axis) - get_coord(r->orig, axis)) * invD;
        if (invD < 0.0) { double tmp = t0; t0 = t1; t1 = tmp; }
        if (t0 > t_min) t_min = t0;
        if (t1 < t_max) t_max = t1;
        if (t_max <= t_min) return 0;
    }
    return 1;
}

struct BVHNode { AABB box; int left; int right; int start; int count; int is_leaf; };

struct Scene {
    SphereLight lights[MAX_LIGHTS]; int num_lights;
    Sphere spheres[MAX_SPHERES]; int num_spheres;
    Triangle triangles[MAX_TRIS]; int num_tris;
    int tri_indices[MAX_TRIS];
    BVHNode bvh_nodes[2 * MAX_TRIS]; int num_bvh_nodes;
};

// --- Basic Hit Functions ---
inline int hit_sphere(const SphereLight* s, int light_index, const Ray* r, double t_min, double t_max, HitRecord* rec) {
    vec3 oc = r->orig - s->center;
    double a = dot(r->dir, r->dir);
    double half_b = dot(oc, r->dir);
    double c = dot(oc, oc) - s->radius * s->radius;
    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return 0;
    double root = (-half_b - sqrt(discriminant)) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrt(discriminant)) / a;
        if (root < t_min || root > t_max) return 0;
    }
    rec->t = root; rec->p = ray_at(r, rec->t);
    vec3 outward = (rec->p - s->center) / s->radius;
    set_face_normal(r, outward, rec);
    rec->mat = &s->mat; rec->emission = s->emission; rec->is_light = 1; rec->light_index = light_index;
    return 1;
}
inline int hit_sphere_geom(const Sphere* s, const Ray* r, double t_min, double t_max, HitRecord* rec) {
    vec3 oc = r->orig - s->center;
    double a = dot(r->dir, r->dir);
    double half_b = dot(oc, r->dir);
    double c = dot(oc, oc) - s->radius * s->radius;
    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return 0;
    double root = (-half_b - sqrt(discriminant)) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrt(discriminant)) / a;
        if (root < t_min || root > t_max) return 0;
    }
    rec->t = root; rec->p = ray_at(r, rec->t);
    vec3 outward = (rec->p - s->center) / s->radius;
    set_face_normal(r, outward, rec);
    rec->mat = &s->mat; rec->emission = vec3(0,0,0); rec->is_light = 0; rec->light_index = -1;
    return 1;
}
inline int hit_triangle(const Triangle* tri, const Ray* r, double t_min, double t_max, HitRecord* rec) {
    const double EPS = 1e-6;
    vec3 v0v1 = tri->v1 - tri->v0; vec3 v0v2 = tri->v2 - tri->v0;
    vec3 pvec = cross(r->dir, v0v2); double det = dot(v0v1, pvec);
    if (fabs(det) < EPS) return 0;
    double invDet = 1.0 / det;
    vec3 tvec = r->orig - tri->v0; double u = dot(tvec, pvec) * invDet;
    if (u < 0.0 || u > 1.0) return 0;
    vec3 qvec = cross(tvec, v0v1); double v = dot(r->dir, qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) return 0;
    double t = dot(v0v2, qvec) * invDet;
    if (t < t_min || t > t_max) return 0;
    rec->t = t; rec->p = ray_at(r, t); rec->normal = tri->normal; 
    rec->front_face = dot(r->dir, tri->normal) < 0; if(!rec->front_face) rec->normal = -rec->normal;
    rec->mat = &tri->mat; rec->emission = vec3(0,0,0); rec->is_light = 0; rec->light_index = -1;
    return 1;
}

// --- BVH Builder (Copied from Original for correctness, SAH optional but keeping it simple/robust) ---
static inline int bvh_build_node(Scene* scene, int start, int end) {
    int node_index = scene->num_bvh_nodes++;
    BVHNode* node = &scene->bvh_nodes[node_index];
    node->start = start; node->count = end - start; node->left = -1; node->right = -1; node->is_leaf = 0;
    AABB bounds = make_empty_box(); AABB centroid_bounds = make_empty_box();
    for (int i = start; i < end; ++i) {
        int tri_idx = scene->tri_indices[i];
        const Triangle* tri = &scene->triangles[tri_idx];
        bounds = surrounding_box(bounds, triangle_bounds(tri));
        centroid_bounds = expand_to_include(centroid_bounds, (tri->v0+tri->v1+tri->v2)/3.0);
    }
    node->box = bounds;
    int n = end - start;
    if (n <= 4) { node->is_leaf = 1; return node_index; }
    vec3 diag = centroid_bounds.maximum - centroid_bounds.minimum;
    int axis = 0;
    if (diag.y() > diag.x() && diag.y() >= diag.z()) axis = 1;
    else if (diag.z() > diag.x() && diag.z() > diag.y()) axis = 2;
    double mid = 0.5 * (get_coord(centroid_bounds.minimum, axis) + get_coord(centroid_bounds.maximum, axis));
    int mid_index = start;
    for (int i = start; i < end; ++i) {
        int tri_idx = scene->tri_indices[i];
        const Triangle* tri = &scene->triangles[tri_idx];
        if (get_coord((tri->v0+tri->v1+tri->v2)/3.0, axis) < mid) {
            int tmp = scene->tri_indices[i]; scene->tri_indices[i] = scene->tri_indices[mid_index]; scene->tri_indices[mid_index] = tmp;
            mid_index++;
        }
    }
    if (mid_index == start || mid_index == end) { node->is_leaf = 1; return node_index; }
    node->left  = bvh_build_node(scene, start, mid_index);
    node->right = bvh_build_node(scene, mid_index, end);
    node->is_leaf = 0;
    return node_index;
}
inline void build_scene_bvh_optimized(Scene* scene) {
    scene->num_bvh_nodes = 0; if (scene->num_tris <= 0) return;
    for (int i = 0; i < scene->num_tris; ++i) scene->tri_indices[i] = i;
    bvh_build_node(scene, 0, scene->num_tris);
}

// --- OPTIMIZATION 1: Shadow Cache & AnyHit ---
// Used for scene_occluded (Shadow Rays)
// 1. Checks last cached blocker first (Temporal coherence optimization)
// 2. Returns immediately on first hit (AnyHit)

// Simple global cache for single-threaded environment
static const void* last_shadow_blocker_sphere = nullptr;
static int last_shadow_blocker_tri_idx = -1;

inline int hit_scene_any(const Scene* scene, const Ray* r, double t_min, double t_max) {
    HitRecord temp;
    
    // 1. Check Cache (Sphere)
    if (last_shadow_blocker_sphere) {
        const Sphere* s = (const Sphere*)last_shadow_blocker_sphere;
        // Verify it's still valid/in scene? Assumed static scene.
        if (hit_sphere_geom(s, r, t_min, t_max, &temp)) return 1;
    }
    
    // 2. Check Cache (Triangle) - HARD to cache specific triangle inside BVH easily without ID.
    // Skipping triangle cache for now to keep BVH simple.

    // 3. Check Spheres (Linear - AnyHit)
    for (int i=0; i<scene->num_spheres; ++i) {
        if (hit_sphere_geom(&scene->spheres[i], r, t_min, t_max, &temp)) {
            last_shadow_blocker_sphere = &scene->spheres[i]; // Update Cache
            return 1;
        }
    }
    
    // 4. Check Lights
    for (int i=0; i<scene->num_lights; ++i) {
        if (hit_sphere(&scene->lights[i], i, r, t_min, t_max, &temp)) return 1;
    }

    // 5. Check Triangles (BVH - AnyHit Logic inline)
    // We reuse the standard BVH but return early? 
    // Implementing a clean stack-less AnyHit BVH traversal is complex.
    // Let's use standard BVH traversal but with early exit flag?
    // Standard bvh_hit_node returns boolean "hit_anything".
    // We just need to make sure we don't search for "closest".
    // Actually, bvh_hit_node in Speed_up.h finds *closest*.
    // We need a specific AnyHit BVH traverser.
    
    // Simple recursive AnyHit BVH
    struct BVHAny {
        static int traverse(const Scene* s, int node_idx, const Ray* r, double t0, double t1) {
            const BVHNode* node = &s->bvh_nodes[node_idx];
            if (!aabb_hit(node->box, r, t0, t1)) return 0;
            if (node->is_leaf) {
                HitRecord tr;
                for(int i=0; i<node->count; ++i) {
                    if(hit_triangle(&s->triangles[s->tri_indices[node->start+i]], r, t0, t1, &tr)) return 1;
                }
                return 0;
            }
            if (node->left >= 0 && traverse(s, node->left, r, t0, t1)) return 1;
            if (node->right >= 0 && traverse(s, node->right, r, t0, t1)) return 1;
            return 0;
        }
    };
    
    if (scene->num_bvh_nodes > 0) {
        if (BVHAny::traverse(scene, 0, r, t_min, t_max)) return 1;
    } else {
        for (int i=0; i<scene->num_tris; ++i) 
            if (hit_triangle(&scene->triangles[i], r, t_min, t_max, &temp)) return 1;
    }

    return 0;
}

// --- Standard Hit Scene (Optimized) ---
// Just aliases to standard logic, but ensures we use the fast math compiler flags
static inline int bvh_hit_node_std(const Scene* scene, int node_index, const Ray* r, double t_min, double t_max, HitRecord* out_rec, double* closest_t) {
    const BVHNode* node = &scene->bvh_nodes[node_index];
    if (!aabb_hit(node->box, r, t_min, t_max)) return 0;
    int hit = 0;
    if (node->is_leaf) {
        HitRecord temp;
        for (int i = 0; i < node->count; ++i) {
            if (hit_triangle(&scene->triangles[scene->tri_indices[node->start + i]], r, t_min, *closest_t, &temp)) {
                *closest_t = temp.t; *out_rec = temp; hit = 1;
            }
        }
    } else {
        if (node->left >= 0) hit |= bvh_hit_node_std(scene, node->left, r, t_min, *closest_t, out_rec, closest_t);
        if (node->right >= 0) hit |= bvh_hit_node_std(scene, node->right, r, t_min, *closest_t, out_rec, closest_t);
    }
    return hit;
}

inline int hit_scene_optimized(const Scene* scene, const Ray* r, double t_min, double t_max, HitRecord* rec) {
    int hit_anything = 0; double closest = t_max;
    if (scene->num_bvh_nodes > 0) {
        HitRecord temp;
        if (bvh_hit_node_std(scene, 0, r, t_min, closest, &temp, &closest)) { hit_anything = 1; *rec = temp; }
    } else {
        HitRecord temp;
        for (int i=0; i<scene->num_tris; ++i) if(hit_triangle(&scene->triangles[i], r, t_min, closest, &temp)) { hit_anything = 1; closest = temp.t; *rec = temp; }
    }
    HitRecord temp;
    for (int i=0; i<scene->num_spheres; ++i) if(hit_sphere_geom(&scene->spheres[i], r, t_min, closest, &temp)) { hit_anything = 1; closest = temp.t; *rec = temp; }
    for (int i=0; i<scene->num_lights; ++i) if(hit_sphere(&scene->lights[i], i, r, t_min, closest, &temp)) { hit_anything = 1; closest = temp.t; *rec = temp; }
    return hit_anything;
}

#endif // SPEED_UP_OPT_H
"""

# ================= 3. 編譯與執行 =================
def generate_source_files():
    print("🛠 [2/5] Generating Source Files & Headers...")
    with open(os.path.join(OUTPUT_DIR, "Speed_up_opt.h"), "w") as f: f.write(CODE_SPEEDUP_OPT)

    try:
        with open("tracer_light_source.h", 'r') as f: light_src = f.read()
        with open("tracer_BSDF_sampling.h", 'r') as f: bsdf_src = f.read()
        with open("tracer_MIS.h", 'r') as f: mis_src = f.read()
    except:
        print("❌ Error: Missing headers."); return False

    # Inject hit_scene_any into scene_occluded
    lo_src = light_src.replace('#include "Speed_up.h"', '#include "Speed_up_opt.h"')
    lo_src = lo_src.replace("return (double)rand() / (double)RAND_MAX;", "return fast_rand01();")
    lo_src = lo_src.replace("hit_scene(", "hit_scene_optimized(")
    lo_src = lo_src.replace("inline int scene_occluded", "inline int scene_occluded_old") 
    new_occluded = """
inline int scene_occluded(const Scene* scene, const point3& p, const vec3& wi, double max_dist) {
    const double EPS = 1e-4;
    Ray shadow_ray; shadow_ray.orig = p + EPS * wi; shadow_ray.dir = wi;
    if (hit_scene_any(scene, &shadow_ray, EPS, max_dist - EPS)) return 1;
    return 0;
}
""" 
    if "#endif" in lo_src: lo_src = lo_src.rsplit("#endif", 1)[0] + new_occluded + "\n#endif"
    else: lo_src += new_occluded
    with open(os.path.join(OUTPUT_DIR, "tracer_light_source_opt.h"), "w") as f: f.write(lo_src)

    # Basic replacements
    bo_src = bsdf_src.replace('#include "Speed_up.h"', '#include "Speed_up_opt.h"')
    bo_src = bo_src.replace('#include "tracer_light_source.h"', '#include "tracer_light_source_opt.h"')
    bo_src = bo_src.replace("rand01()", "fast_rand01()")
    with open(os.path.join(OUTPUT_DIR, "tracer_BSDF_sampling_opt.h"), "w") as f: f.write(bo_src)

    # Inject optimized includes
    mo_src = mis_src.replace('#include "tracer_BSDF_sampling.h"', '#include "tracer_BSDF_sampling_opt.h"')
    mo_src = mo_src.replace("rand01()", "fast_rand01()")
    mo_src = mo_src.replace("hit_scene(", "hit_scene_optimized(")
    # Keep original ray_color (no RR injection for visual identity check)
    # But we use hit_scene_optimized which is faster.
    with open(os.path.join(OUTPUT_DIR, "tracer_MIS_opt.h"), "w") as f: f.write(mo_src)
    return True

def compile_and_run():
    print("🔨 [3/5] Compiling...")
    if not os.path.exists(SOURCE_FILE): return

    with open(SOURCE_FILE, 'r') as f: src = f.read()
    src = re.sub(r'(chdir\s*\(\s*dirname\s*\(\s*exe_path\s*\)\s*\)\s*;)', r'// \1', src)
    
    # Filter to only MIS
    src_clean = re.sub(r'(render_image\s*\([^;\{]+AdvCG_light\.ppm[^;\{]+\);)', r'/* \1 */', src)
    src_clean = re.sub(r'(render_image\s*\([^;\{]+AdvCG_bsdf\.ppm[^;\{]+\);)', r'/* \1 */', src_clean)
    src_clean = re.sub(r'(render_mis_weight_image\s*\([^;\{]+\);)', r'/* \1 */', src_clean)
    
    # 1. Original
    if os.path.exists(HEADER_FILE):
        with open(HEADER_FILE, 'r') as f: orig_h = f.read()
        with open(os.path.join(OUTPUT_DIR, "Speed_up.h"), 'w') as f: f.write(orig_h)
    with open(os.path.join(OUTPUT_DIR, "Original.cc"), 'w') as f: f.write(src_clean)
    
    # 2. Optimized
    src_opt = src_clean
    src_opt = src_opt.replace('#include "Speed_up.h"', '#include "Speed_up_opt.h"')
    src_opt = src_opt.replace('#include "tracer_MIS.h"', '#include "tracer_MIS_opt.h"')
    src_opt = src_opt.replace('build_scene_bvh(&scene);', 'build_scene_bvh_optimized(&scene);')
    src_opt = src_opt.replace('build_scene_bvh(scene);', 'build_scene_bvh_optimized(scene);')
    src_opt = src_opt.replace('rand01()', 'fast_rand01()')
    with open(os.path.join(OUTPUT_DIR, "Optimized.cc"), 'w') as f: f.write(src_opt)

    flags_orig = ["-O3", "-std=c++11", f"-I{OUTPUT_DIR}", "-I.", '-DOUT_DIR=""']
    flags_opt  = ["-O3", "-std=c++11", "-ffast-math", f"-I{OUTPUT_DIR}", "-I.", '-DOUT_DIR=""']

    print("   [Compile] Original...")
    subprocess.run([COMPILER] + flags_orig + [os.path.join(OUTPUT_DIR, "Original.cc"), "-o", EXE_FAST], check=True)
    
    print("   [Compile] Optimized...")
    subprocess.run([COMPILER] + flags_opt + [os.path.join(OUTPUT_DIR, "Optimized.cc"), "-o", EXE_OPT], check=True)

    print("\n🏁 [4/5] Running Benchmark (1024x1024, Only MIS)...")
    
    print("\n>>> Running Original (Before Optimization)...")
    t0 = time.time()
    subprocess.run([os.path.abspath(EXE_FAST), "input_optimized.txt"], cwd=OUTPUT_DIR)
    t_orig = time.time() - t0
    if os.path.exists(os.path.join(OUTPUT_DIR, "AdvCG_mis.ppm")):
        shutil.move(os.path.join(OUTPUT_DIR, "AdvCG_mis.ppm"), os.path.join(OUTPUT_DIR, "before_optimization.ppm"))
    print(f"⏱️  Original Time: {t_orig:.4f}s")

    print("\n>>> Running Optimized (After Optimization)...")
    t0 = time.time()
    subprocess.run([os.path.abspath(EXE_OPT), "input_optimized.txt"], cwd=OUTPUT_DIR)
    t_opt = time.time() - t0
    if os.path.exists(os.path.join(OUTPUT_DIR, "AdvCG_mis.ppm")):
        shutil.move(os.path.join(OUTPUT_DIR, "AdvCG_mis.ppm"), os.path.join(OUTPUT_DIR, "after_optimization.ppm"))
    print(f"⏱️  Optimized Time: {t_opt:.4f}s")

    print("\n=== [V17 Results] ===")
    print(f"Original : {t_orig:.4f}s")
    print(f"Optimized: {t_opt:.4f}s")
    speedup = t_orig / t_opt if t_opt > 0 else 0
    print(f"🚀 Speedup: {speedup:.2f}x")
    print(f"🖼  Images Generated: before_optimization.ppm, after_optimization.ppm")

if __name__ == "__main__":
    if prepare_scene():
        if generate_source_files():
            compile_and_run()