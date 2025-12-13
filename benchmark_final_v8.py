import os
import time
import subprocess
import re
import math
import platform
import sys

# ================= 參數設定 =================
OUTPUT_DIR = "benchmark_output"
COMPILER = "g++"
SOURCE_FILE = "AdvCG_Final_MIS.cc"
HEADER_FILE = "Speed_up.h"

# 來源檔案 (你的原始場景)
SOURCE_INPUT = "input.txt"
# 目標檔案 (轉換後的高密度網格場景)
TARGET_INPUT = os.path.join(OUTPUT_DIR, "input_mesh.txt")

# 設定球體細分度
# 50 steps => 每顆球約 5000 面
# 你的場景有 ~200 顆球 => 總共 ~100 萬面 => 暴力法會非常慢
SPHERE_STEPS = 50 

if len(sys.argv) > 1 and sys.argv[1] == "quick":
    print("⚡️ [Quick Mode] 啟動快速測試 (低細分)...")
    SPHERE_STEPS = 10 # 快速測試用
else:
    print("🐢 [Full Mode] 啟動完整壓力測試 (高細分)...")

if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
exe_ext = ".exe" if platform.system() == "Windows" else ""
EXE_FAST = os.path.join(OUTPUT_DIR, f"render_fast{exe_ext}")
EXE_SLOW = os.path.join(OUTPUT_DIR, f"render_slow{exe_ext}")

# ================= 1. 強制修補 Header =================
def fix_header():
    print(f"🔧 [1/4] 強制修補 Speed_up.h (MAX_TRIS=2000000)...")
    if not os.path.exists(HEADER_FILE):
        print(f"❌ 找不到 {HEADER_FILE}")
        return False
    with open(HEADER_FILE, 'r', encoding='utf-8') as f: content = f.read()
    # 開 200 萬比較保險
    new_content = re.sub(r'#define\s+MAX_TRIS\s+\d+', f'#define MAX_TRIS 2000000', content)
    with open(os.path.join(OUTPUT_DIR, HEADER_FILE), 'w', encoding='utf-8') as f: f.write(new_content)
    return True

# ================= 2. 轉換場景 (input.txt -> input_mesh.txt) =================
def convert_scene():
    if not os.path.exists(SOURCE_INPUT):
        print(f"❌ 錯誤：找不到 {SOURCE_INPUT}，請確認檔案在同一目錄下。")
        return False

    print(f"⚙️ [2/4] 正在轉換場景 {SOURCE_INPUT} -> {TARGET_INPUT} ...")
    
    total_tris = 0
    sphere_count = 0

    def get_sphere_mesh(cx, cy, cz, r, steps):
        tris = []
        for i in range(steps):
            lat0, lat1 = math.pi * (-0.5 + i/steps), math.pi * (-0.5 + (i+1)/steps)
            z0, zr0 = math.sin(lat0), math.cos(lat0)
            z1, zr1 = math.sin(lat1), math.cos(lat1)
            for j in range(steps):
                lng0, lng1 = 2*math.pi * j/steps, 2*math.pi * (j+1)/steps
                x0, y0 = math.cos(lng0), math.sin(lng0)
                x1, y1 = math.cos(lng1), math.sin(lng1)
                p1 = (cx+r*x0*zr0, cy+r*y0*zr0, cz+r*z0)
                p2 = (cx+r*x0*zr1, cy+r*y0*zr1, cz+r*z1)
                p3 = (cx+r*x1*zr0, cy+r*y1*zr0, cz+r*z0)
                p4 = (cx+r*x1*zr1, cy+r*y1*zr1, cz+r*z1)
                tris.append((p1,p2,p3)); tris.append((p3,p2,p4))
        return tris

    with open(SOURCE_INPUT, 'r', encoding='utf-8') as f_in, \
         open(TARGET_INPUT, 'w', encoding='utf-8') as f_out:
        
        for line in f_in:
            parts = line.split()
            if not parts:
                continue
            
            tag = parts[0]
            
            # 遇到球體 S，轉換成 T
            if tag == 'S':
                try:
                    x, y, z, r = map(float, parts[1:5])
                    mesh = get_sphere_mesh(x, y, z, r, SPHERE_STEPS)
                    f_out.write(f"# --- Converted Sphere {sphere_count} ({len(mesh)} tris) ---\n")
                    for t in mesh:
                        p1, p2, p3 = t
                        f_out.write(f"T {p1[0]:.4f} {p1[1]:.4f} {p1[2]:.4f} {p2[0]:.4f} {p2[1]:.4f} {p2[2]:.4f} {p3[0]:.4f} {p3[1]:.4f} {p3[2]:.4f}\n")
                        total_tris += 1
                    sphere_count += 1
                except ValueError:
                    f_out.write(line) # 解析失敗保留原樣
            else:
                # 其他 E, V, M, SL 等保留原樣
                f_out.write(line)

    print(f"✅ 轉換完成：共 {sphere_count} 顆球 -> {total_tris} 個三角形。")
    return True

# ================= 3. 編譯與執行 =================
def compile_and_run():
    print("🔨 [3/4] 處理原始碼 (只保留 MIS Render)...")
    with open(SOURCE_FILE, 'r') as f: src = f.read()
    
    # 修正 chdir & stack
    src = re.sub(r'(chdir\s*\(\s*dirname\s*\(\s*exe_path\s*\)\s*\)\s*;)', r'// \1', src)
    src = re.sub(r'(\s+)Scene\s+scene;', r'\1// Scene scene;\n\1static Scene scene;', src)
    
    # 插入診斷
    src = re.sub(r'(if\s*\(!load_scene_from_file[^\{]+\{)', r'\1', src)
    src = src.replace('fprintf(stderr, "Scene loading failed.\\n");\n        return 1;\n    }', 'fprintf(stderr, "Scene loading failed.\\n"); return 1; } printf("[DEBUG] Tris loaded: %d\\n", scene.num_tris);')

    # 【關鍵修正】只保留 AdvCG_mis.ppm
    # 使用 /* */ 註解掉其他 render 呼叫
    src = re.sub(r'(render_image\s*\([^;\{]+AdvCG_light\.ppm[^;\{]+\);)', r'/* \1 */', src)
    src = re.sub(r'(render_image\s*\([^;\{]+AdvCG_bsdf\.ppm[^;\{]+\);)', r'/* \1 */', src)
    src = re.sub(r'(render_mis_weight_image\s*\([^;\{]+\);)', r'/* \1 */', src)

    src_fast = src
    src_slow = re.sub(r'(build_scene_bvh\s*\(\s*scene\s*\)\s*;)', r'// \1', src)

    p_fast = os.path.join(OUTPUT_DIR, "Fast.cc")
    p_slow = os.path.join(OUTPUT_DIR, "Slow.cc")
    with open(p_fast, 'w') as f: f.write(src_fast)
    with open(p_slow, 'w') as f: f.write(src_slow)

    flags = ["-O3", "-std=c++11", f"-I{OUTPUT_DIR}", "-I.", '-DOUT_DIR=""']
    
    print("   Compiling...")
    subprocess.run([COMPILER] + flags + [p_fast, "-o", EXE_FAST], check=True)
    subprocess.run([COMPILER] + flags + [p_slow, "-o", EXE_SLOW], check=True)

    print("\n🚀 [4/4] 開始跑分...")
    abs_exe_slow = os.path.abspath(EXE_SLOW)
    abs_exe_fast = os.path.abspath(EXE_FAST)
    
    print("\n>>> Running No BVH (Linear)...")
    t0 = time.time()
    subprocess.run([abs_exe_slow, os.path.basename(TARGET_INPUT)], cwd=OUTPUT_DIR, check=True)
    t_slow = time.time() - t0
    print(f">>> Time: {t_slow:.4f}s")

    print("\n>>> Running With BVH (Tree)...")
    t0 = time.time()
    subprocess.run([abs_exe_fast, os.path.basename(TARGET_INPUT)], cwd=OUTPUT_DIR, check=True)
    t_fast = time.time() - t0
    print(f">>> Time: {t_fast:.4f}s")

    print("\n=== [結果] ===")
    print(f"No BVH: {t_slow:.4f}s")
    print(f"BVH   : {t_fast:.4f}s")
    if t_fast > 0:
        print(f"🚀 Speedup: {t_slow/t_fast:.2f}x")
        
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6,4))
        plt.bar(['No BVH', 'With BVH'], [t_slow, t_fast], color=['#ff9999', '#66b3ff'])
        plt.title(f'Speedup: {t_slow/t_fast:.1f}x (Mesh converted from input.txt)')
        plt.ylabel('Time (s)')
        for i, v in enumerate([t_slow, t_fast]):
            plt.text(i, v, f"{v:.2f}s", ha='center', va='bottom')
        plt.savefig(os.path.join(OUTPUT_DIR, "benchmark_chart.png"))
        print(f"📊 圖表已存: {os.path.join(OUTPUT_DIR, 'benchmark_chart.png')}")
    except: pass

if __name__ == "__main__":
    if fix_header():
        if convert_scene():
            compile_and_run()