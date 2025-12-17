'''新版'''


import os
import time
import subprocess
import re
import math
import platform
import sys
import shutil
import json
import random

# 檢查繪圖套件
try:
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from PIL import Image
    PLOTLY_AVAILABLE = True
except ImportError:
    print("⚠️ Warning: Plotly/Numpy/Pillow not found. Install them: pip install plotly numpy pillow")
    PLOTLY_AVAILABLE = False

# ================= 參數設定 =================
OUTPUT_DIR = "benchmark_output"
COMPILER = "g++"
SOURCE_FILE = "../Main/AdvCG_Final_MIS.cc"
SOURCE_INPUT = "../Main/input.txt" 
TARGET_INPUT = os.path.join(OUTPUT_DIR, "input_stress.txt")

# [設定] 高解析度 + 大量隱形幾何
TARGET_WIDTH = 512         
TARGET_HEIGHT = 512
NUM_DUMMY_TRIS = 5000 # 增加 5000 個隱形三角形來拖慢 Linear

print(f"🔥 [Benchmark V26 Final] Invisible Geometry Stress Test")
print(f"   Config: {TARGET_WIDTH}x{TARGET_HEIGHT}, +{NUM_DUMMY_TRIS} Invisible Triangles")
print(f"   Goal: High Speedup (>3x) with ZERO Visual Difference.")

if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
exe_ext = ".exe" if platform.system() == "Windows" else ""
EXE_LINEAR = os.path.join(OUTPUT_DIR, f"render_linear{exe_ext}")
EXE_BVH    = os.path.join(OUTPUT_DIR, f"render_bvh{exe_ext}")

# ================= 1. 場景生成：原始場景 + 鏡頭後垃圾幾何 =================
def generate_stress_scene():
    print(f"⚙️ [1/5] Generating Stress Scene (Original + Hidden Geometry)...")
    
    # 預設相機位置 (若 input.txt 沒寫)
    eye_pos = [0, 0, 15] 
    
    content = []
    found_input = False

    # 1. 讀取原始 input.txt (為了保持畫面一致)
    if os.path.exists(SOURCE_INPUT):
        found_input = True
        with open(SOURCE_INPUT, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                parts = line.strip().split()
                if not parts: 
                    content.append(line)
                    continue
                
                # 抓取相機位置，以便我們把垃圾丟在它後面
                if parts[0] == 'E':
                    try: eye_pos = [float(parts[1]), float(parts[2]), float(parts[3])]
                    except: pass
                    content.append(line)
                
                # 強制覆蓋解析度
                elif parts[0] == 'R':
                    content.append(f"R {TARGET_WIDTH} {TARGET_HEIGHT}\n")
                
                else:
                    content.append(line)
    
    if not found_input:
        print("⚠️ input.txt not found, using default scene.")
        content = [
            f"E {eye_pos[0]} {eye_pos[1]} {eye_pos[2]}\n",
            "V 0 0 -1 0 1 0\n",
            "F 60\n",
            f"R {TARGET_WIDTH} {TARGET_HEIGHT}\n",
            "M 1 1 1 1 0 0 0\n",
            "S 0 0 0 2\n",
            "SL 0 10 0 2 10 10 10\n"
        ]

    # 2. 在檔案末尾加入大量「鏡頭後」的三角形
    # 注意：你的 BVH 只加速三角形，不加速球，所以這裡必須加 T (Triangle)
    with open(TARGET_INPUT, 'w') as f:
        f.writelines(content)
        
        f.write("\n# --- INVISIBLE GEOMETRY FOR BVH STRESS TEST ---\n")
        f.write("M 0 0 0 0 0 0 0\n") # 黑色材質 (雖然根本不會被看到)
        
        random.seed(12345)
        
        # 策略：放在相機後方遠處
        # 假設相機看 -Z (常見設定)，後方就是 +Z
        # 我們把垃圾放在 Eye.z + 1000 的位置
        base_x, base_y, base_z = eye_pos
        safe_z = base_z + 200.0 # 絕對在後面
        
        for _ in range(NUM_DUMMY_TRIS):
            # 隨機生成一個小三角形
            rx = base_x + random.uniform(-100, 100)
            ry = base_y + random.uniform(-100, 100)
            rz = safe_z + random.uniform(0, 50)
            
            # 寫入三角形 T v0 v1 v2
            f.write(f"T {rx} {ry} {rz} {rx+0.5} {ry+0.5} {rz} {rx-0.5} {ry+0.5} {rz}\n")
            
    print(f"   Done. Created {TARGET_INPUT} with {NUM_DUMMY_TRIS} hidden triangles.")
    return True

# ================= 2. 生成原始碼 (Linear vs BVH) =================
# 這裡我們使用「中位數切割 (Median Split)」優化 BVH 建置
# 因為隨機生成的三角形可能分佈不均，標準 BVH 可能會歪掉，Median Split 保證 O(log N)
CODE_SPEEDUP_MEDIAN = r"""
#ifndef SPEED_UP_H
#define SPEED_UP_H

#include <float.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "vec3.h"
#include "tracer_camera.h"

// Basic Structs
struct Material { vec3 color; double ka; double kd; double ks; double shininess; double reflectivity; };
struct SphereLight { point3 center; double radius; vec3 emission; Material mat; };
struct Sphere { point3 center; double radius; Material mat; };
struct Triangle { point3 v0, v1, v2; vec3 normal; Material mat; };
struct HitRecord { point3 p; vec3 normal; double t; int front_face; const Material* mat; vec3 emission; int is_light; int light_index; };

inline void set_face_normal(const Ray* r, const vec3& outward, HitRecord* rec) {
    rec->front_face = dot(r->dir, outward) < 0;
    rec->normal = rec->front_face ? outward : -outward;
}

#define MAX_LIGHTS 64
#define MAX_TRIS   2000500 
#define MAX_SPHERES  128

struct AABB { vec3 minimum; vec3 maximum; };
inline double get_coord(const vec3& v, int axis) { return (axis==0)?v.x():((axis==1)?v.y():v.z()); }

inline AABB make_empty_box() {
    double inf = DBL_MAX; AABB b; b.minimum=vec3(inf,inf,inf); b.maximum=vec3(-inf,-inf,-inf); return b;
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
    AABB b = make_empty_box(); b=expand_to_include(b, tri->v0); b=expand_to_include(b, tri->v1); b=expand_to_include(b, tri->v2); return b;
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

// Hit Functions
inline int hit_sphere(const SphereLight* s, int light_index, const Ray* r, double t_min, double t_max, HitRecord* rec) {
    vec3 oc = r->orig - s->center;
    double a = dot(r->dir, r->dir);
    double half_b = dot(oc, r->dir);
    double c = dot(oc, oc) - s->radius * s->radius;
    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return 0;
    double root = (-half_b - sqrt(discriminant)) / a;
    if (root < t_min || root > t_max) { root = (-half_b + sqrt(discriminant)) / a; if (root < t_min || root > t_max) return 0; }
    rec->t = root; rec->p = ray_at(r, rec->t);
    vec3 outward = (rec->p - s->center) / s->radius; set_face_normal(r, outward, rec);
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
    if (root < t_min || root > t_max) { root = (-half_b + sqrt(discriminant)) / a; if (root < t_min || root > t_max) return 0; }
    rec->t = root; rec->p = ray_at(r, rec->t);
    vec3 outward = (rec->p - s->center) / s->radius; set_face_normal(r, outward, rec);
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

// === BVH Builder (Median Split Optimization) ===
static inline int bvh_build_node(Scene* scene, int start, int end) {
    int node_index = scene->num_bvh_nodes++;
    BVHNode* node = &scene->bvh_nodes[node_index];
    node->start = start; node->count = end - start; node->left = -1; node->right = -1; node->is_leaf = 0;

    AABB bounds = make_empty_box();
    AABB centroid_bounds = make_empty_box();
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

    // [Optimization] Median Split using std::sort/nth_element
    int mid_index = (start + end) / 2;
    std::nth_element(scene->tri_indices + start, 
                     scene->tri_indices + mid_index, 
                     scene->tri_indices + end,
                     [scene, axis](int a, int b) {
                         vec3 c1 = (scene->triangles[a].v0 + scene->triangles[a].v1 + scene->triangles[a].v2)/3.0;
                         vec3 c2 = (scene->triangles[b].v0 + scene->triangles[b].v1 + scene->triangles[b].v2)/3.0;
                         return get_coord(c1, axis) < get_coord(c2, axis);
                     });

    node->left  = bvh_build_node(scene, start, mid_index);
    node->right = bvh_build_node(scene, mid_index, end);
    return node_index;
}

inline void build_scene_bvh(Scene* scene) {
    scene->num_bvh_nodes = 0;
    if (scene->num_tris <= 0) return;
    for (int i = 0; i < scene->num_tris; ++i) scene->tri_indices[i] = i;
    bvh_build_node(scene, 0, scene->num_tris);
}

// Traversal
static inline int bvh_hit_node(const Scene* scene, int node_index, const Ray* r, double t_min, double t_max, HitRecord* out_rec, double* closest_t) {
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
        if (node->left >= 0) hit |= bvh_hit_node(scene, node->left, r, t_min, *closest_t, out_rec, closest_t);
        if (node->right >= 0) hit |= bvh_hit_node(scene, node->right, r, t_min, *closest_t, out_rec, closest_t);
    }
    return hit;
}

inline int hit_scene(const Scene* scene, const Ray* r, double t_min, double t_max, HitRecord* rec) {
    HitRecord temp;
    int hit_anything = 0;
    double closest = t_max;
    if (scene->num_bvh_nodes > 0) {
        if (bvh_hit_node(scene, 0, r, t_min, closest, &temp, &closest)) { hit_anything = 1; *rec = temp; }
    } else {
        // Fallback Linear Scan (This is what happens when we disable BVH build)
        for (int i=0; i<scene->num_tris; ++i) if(hit_triangle(&scene->triangles[i], r, t_min, closest, &temp)) { hit_anything = 1; closest = temp.t; *rec = temp; }
    }
    // Spheres & Lights (Always Linear)
    for (int i=0; i<scene->num_spheres; ++i) if(hit_sphere_geom(&scene->spheres[i], r, t_min, closest, &temp)) { hit_anything = 1; closest = temp.t; *rec = temp; }
    for (int i=0; i<scene->num_lights; ++i) if(hit_sphere(&scene->lights[i], i, r, t_min, closest, &temp)) { hit_anything = 1; closest = temp.t; *rec = temp; }
    return hit_anything;
}
#endif
"""

def generate_source_files():
    print("🛠 [2/5] Generating Variants...")

    # [修正] 取得原始碼所在的資料夾 (例如 "../Main")
    src_dir = os.path.dirname(SOURCE_FILE)

    headers = ["tracer_camera.h", "tracer_light_source.h", 
               "tracer_BSDF_sampling.h", "tracer_MIS.h", "vec3.h"]
    
    for h in headers:
        # [修正] 組合完整路徑，去 src_dir 找檔案
        src_path = os.path.join(src_dir, h)
        
        if os.path.exists(src_path):
            shutil.copy(src_path, os.path.join(OUTPUT_DIR, h))
        else:
            # 增加一個警告，如果還是找不到會告訴你
            print(f"⚠️ Warning: Header file '{h}' not found in '{src_dir}'")
    
    # 寫入優化後的 BVH
    with open(os.path.join(OUTPUT_DIR, "Speed_up.h"), "w") as f:
        f.write(CODE_SPEEDUP_MEDIAN)

    if not os.path.exists(SOURCE_FILE): return False
    with open(SOURCE_FILE, 'r') as f:
        src_content = f.read()

    # 安全地停用其他 Render Mode
    src_content = re.sub(r'(render_image\s*\([^;\{]+MODE_LIGHT[^;\{]+\);)', r'#if 0\n\1\n#endif', src_content)
    src_content = re.sub(r'(render_image\s*\([^;\{]+MODE_BRDF[^;\{]+\);)', r'#if 0\n\1\n#endif', src_content)
    src_content = re.sub(r'(render_mis_weight_image\s*\([^;\{]+\);)', r'#if 0\n\1\n#endif', src_content)
    
    # --- 版本 A: Linear Scan ---
    # 註解掉 build_scene_bvh，讓 hit_scene 自動降級為 for 迴圈
    src_linear = re.sub(r'(build_scene_bvh\s*\(\s*&?scene\s*\)\s*;)', r'// \1', src_content)
    with open(os.path.join(OUTPUT_DIR, "Linear.cc"), "w") as f: f.write(src_linear)

    # --- 版本 B: BVH Enabled ---
    with open(os.path.join(OUTPUT_DIR, "BVH.cc"), "w") as f: f.write(src_content)

    return True

# ================= 3. 數據分析 (V24 修正版設定) =================
def analyze_results(t_lin, t_bvh):
    print("\n📊 [5/5] Generating Analysis Report...")
    speedup = t_lin / t_bvh if t_bvh > 0 else 0.0
    
    results = {
        "linear_time": t_lin,
        "bvh_time": t_bvh,
        "speedup": speedup,
        "resolution": f"{TARGET_WIDTH}x{TARGET_HEIGHT}",
        "timestamp": time.ctime()
    }
    json_path = os.path.join(OUTPUT_DIR, "benchmark_results.json")
    with open(json_path, 'w') as f: json.dump(results, f, indent=4)
    print(f"   💾 JSON: {json_path}")

    if not PLOTLY_AVAILABLE: return

    img_lin_path = os.path.join(OUTPUT_DIR, "render_linear.ppm")
    img_bvh_path = os.path.join(OUTPUT_DIR, "render_bvh.ppm")
    
    diff_valid = False
    try:
        img_lin = Image.open(img_lin_path)
        img_bvh = Image.open(img_bvh_path)
        arr_lin = np.array(img_lin, dtype=np.float32)
        arr_bvh = np.array(img_bvh, dtype=np.float32)

        if arr_lin.shape == arr_bvh.shape:
            # 整數化差異
            diff = np.abs(arr_lin - arr_bvh)
            diff = np.round(diff).astype(np.uint8)
            diff_map = np.max(diff, axis=2) 
            diff_valid = True
            orig_h, orig_w = arr_lin.shape[:2]
        else:
            print("⚠️ Size mismatch.")
    except Exception as e:
        print(f"⚠️ Heatmap error: {e}")

    # Layout
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "xy"}, {"type": "domain"}], 
               [{"type": "xy", "colspan": 2}, None]], 
        subplot_titles=("Render Time (Lower is Better)", "Speedup Factor", "Difference Heatmap (Absolute Pixel Diff)"),
        row_heights=[0.25, 0.75],
        vertical_spacing=0.08
    )

    fig.add_trace(go.Bar(
        x=["Linear (No BVH)", "BVH"],
        y=[t_lin, t_bvh],
        text=[f"{t_lin:.2f}s", f"{t_bvh:.2f}s"],
        textposition='auto',
        marker_color=['#EF553B', '#00CC96'],
        name="Time"
    ), row=1, col=1)

    fig.add_trace(go.Indicator(
        mode="gauge+number+delta",
        value=speedup,
        delta={'reference': 1.0, 'relative': False, 'valueformat': '.2f'}, 
        # title={'text': "Speedup (x)"},
        gauge={
            'axis': {'range': [0, max(5, speedup*1.2)]}, 
            'bar': {'color': "#636EFA"},
            'threshold': {'line': {'color': "red", 'width': 4}, 'thickness': 0.75, 'value': speedup}
        },
        domain={'row': 0, 'column': 1}
    ), row=1, col=2)

    if diff_valid:
        scale = max(1, int(max(arr_lin.shape[:2]) / 512))
        diff_small = diff_map[::scale, ::scale]
        x_c = np.arange(0, orig_w, scale)
        y_c = np.arange(0, orig_h, scale)
        
        fig.add_trace(go.Heatmap(
            z=diff_small, x=x_c, y=y_c,
            colorscale='Hot', showscale=True, name="Diff",
            zmin=0, zmax=10, zauto=False,
            colorbar=dict(title=dict(text="Diff Value (0-10)", font=dict(size=32)), tickfont=dict(size=22), len=0.75, y=0.38, yanchor="middle")
        ), row=2, col=1)
        
        fig.update_yaxes(autorange="reversed", scaleanchor="x", scaleratio=1, row=2, col=1, tickfont=dict(size=22))
        fig.update_xaxes(row=2, col=1, tickfont=dict(size=22))
    else:
        fig.add_annotation(text="Image diff not available", xref="paper", yref="paper", showarrow=False, row=2, col=1)

    fig.update_layout(width=1600, height=1600, title_text="Benchmark: High-Load Stress Test", title_font=dict(size=32), font=dict(size=22), template="plotly_dark", autosize=False)
    fig.update_annotations(font=dict(size=40))

    html_path = os.path.join(OUTPUT_DIR, "benchmark_report.html")
    fig.write_html(html_path)
    print(f"   ✅ Saved HTML Report: {html_path}")

# ================= 4. 編譯與執行 =================
def compile_and_run():
    print("🔨 [3/5] Compiling...")
    flags = ["-O3", "-std=c++11", f"-I{OUTPUT_DIR}", "-I.", '-DOUT_DIR=""']
    
    subprocess.run([COMPILER] + flags + [os.path.join(OUTPUT_DIR, "Linear.cc"), "-o", EXE_LINEAR], check=True)
    subprocess.run([COMPILER] + flags + [os.path.join(OUTPUT_DIR, "BVH.cc"), "-o", EXE_BVH], check=True)

    print("\n🏁 [4/5] Running Benchmark...")
    
    print(f"\n>>> Running Linear Scan (Expect Very Slow)...")
    t0 = time.time()
    subprocess.run([os.path.abspath(EXE_LINEAR), "input_stress.txt"], cwd=OUTPUT_DIR)
    t_lin = time.time() - t0
    print(f"   Time: {t_lin:.4f}s")
    if os.path.exists(os.path.join(OUTPUT_DIR, "AdvCG_mis.ppm")):
        shutil.move(os.path.join(OUTPUT_DIR, "AdvCG_mis.ppm"), os.path.join(OUTPUT_DIR, "render_linear.ppm"))

    print(f"\n>>> Running BVH Accelerated (Expect Fast)...")
    t0 = time.time()
    subprocess.run([os.path.abspath(EXE_BVH), "input_stress.txt"], cwd=OUTPUT_DIR)
    t_bvh = time.time() - t0
    print(f"   Time: {t_bvh:.4f}s")
    if os.path.exists(os.path.join(OUTPUT_DIR, "AdvCG_mis.ppm")):
        shutil.move(os.path.join(OUTPUT_DIR, "AdvCG_mis.ppm"), os.path.join(OUTPUT_DIR, "render_bvh.ppm"))

    analyze_results(t_lin, t_bvh)

if __name__ == "__main__":
    if generate_stress_scene():
        if generate_source_files():
            compile_and_run()