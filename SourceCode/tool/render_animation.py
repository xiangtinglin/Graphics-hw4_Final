''' 
video:
ffmpeg -framerate 24 -i animation_frames_original/frame_%03d.ppm -c:v libx264 -pix_fmt yuv420p ./animation_frames_original/output.mp4
'''
 
import os
import subprocess
import math
import shutil
import sys
import re
import time

# ================= 動畫設定 =================
TOTAL_FRAMES = 24        # 24 幀 (約1秒)
RESOLUTION = 512         # 解析度
ANIM_SPEED = 2.0         # 動作速度
LIGHT_DESCENT_SPEED = 0.2 # 光源移動速度

# 路徑設定
OUTPUT_DIR = "animation_frames_original"  # 改個名字避免跟優化版混淆
BUILD_DIR = "animation_build_original"
COMPILER = "g++"
SOURCE_FILE = "AdvCG_Final_MIS.cc"
INPUT_FILE = "input.txt"

# 確保當前目錄正確
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
if not os.path.exists(BUILD_DIR): os.makedirs(BUILD_DIR)

EXE_PATH = os.path.join(BUILD_DIR, "renderer_anim_original")
if sys.platform == "win32": EXE_PATH += ".exe"

# ================= 1. 編譯原始渲染器 (無優化) =================
def compile_renderer():
    print("🔨 Compiling ORIGINAL renderer (No Optimization)...")
    if not os.path.exists(SOURCE_FILE):
        print(f"❌ Source file {SOURCE_FILE} not found.")
        return False

    with open(SOURCE_FILE, 'r') as f: src = f.read()
    
    # 為了自動化，我們還是需要修改一下原始碼來只渲染 MIS，避免浪費時間
    # 除此之外，不改動任何演算法邏輯
    
    src = re.sub(r'(render_image\s*\([^;\{]+AdvCG_light\.ppm[^;\{]+\);)', r'/* \1 */', src)
    src = re.sub(r'(render_image\s*\([^;\{]+AdvCG_bsdf\.ppm[^;\{]+\);)', r'/* \1 */', src)
    src = re.sub(r'(render_mis_weight_image\s*\([^;\{]+\);)', r'/* \1 */', src)

    # 確保使用原始的 Speed_up.h，不替換
    # (預設 AdvCG_Final_MIS.cc 就是 include "Speed_up.h")

    src_path = os.path.join(BUILD_DIR, "anim_main_original.cc")
    with open(src_path, 'w') as f: f.write(src)

    # 複製原始 header 到 build 目錄，確保讀到正確的檔案
    headers = ["Speed_up.h", "tracer_light_source.h", "tracer_BSDF_sampling.h", "tracer_MIS.h", "vec3.h", "tracer_camera.h"]
    for h in headers:
        if os.path.exists(h):
            shutil.copy(h, os.path.join(BUILD_DIR, h))
        else:
            print(f"❌ Missing header: {h}")
            return False

    # 編譯參數 (只開基本優化 -O3，不開 fast-math 或其他特殊 flag)
    flags = ["-O3", "-std=c++11", f"-I{BUILD_DIR}", f"-I{SCRIPT_DIR}", '-DOUT_DIR=""']
    
    try:
        subprocess.check_call([COMPILER] + flags + ["-o", EXE_PATH, src_path])
        print("✅ Compilation successful!")
    except subprocess.CalledProcessError:
        print("❌ Compilation Failed.")
        return False
    return True

# ================= 2. 產生幀 (與之前邏輯相同) =================
def generate_frame_input(frame_idx):
    if not os.path.exists(INPUT_FILE):
        print(f"❌ Input file {INPUT_FILE} not found!")
        return None

    # 讀取原始 input.txt
    with open(INPUT_FILE, 'r') as f:
        lines = f.readlines()

    output_lines = []
    
    for line in lines:
        parts = line.split()
        if not parts:
            output_lines.append(line)
            continue
            
        # 1. 攔截解析度
        if parts[0] == 'R':
            output_lines.append(f"R {RESOLUTION} {RESOLUTION}\n")
            
        # 2. 光源 (SL) -> Z 軸遞減
        elif parts[0] == 'SL':
            try:
                x, y, z, r = map(float, parts[1:5])
                new_z = z + (frame_idx * LIGHT_DESCENT_SPEED)
                rest_of_line = " ".join(parts[5:])
                output_lines.append(f"SL {x:.3f} {y:.3f} {new_z:.3f} {r:.3f} {rest_of_line}\n")
            except:
                output_lines.append(line)
                
        # 3. 球體 (S) -> Y 軸浮動 (這是玻璃球)
        elif parts[0] == 'S':
            try:
                x, y, z, r = map(float, parts[1:5])
                # 簡單的上下浮動
                offset = math.sin(frame_idx * 0.5 + x) * 0.5 
                new_y = y + offset
                
                rest_of_line = ""
                if len(parts) > 5:
                    rest_of_line = " " + " ".join(parts[5:])
                output_lines.append(f"S {x:.3f} {new_y:.3f} {z:.3f} {r:.3f}{rest_of_line}\n")
            except:
                output_lines.append(line)
        else:
            output_lines.append(line)

    # 寫入 input 檔
    filename = os.path.join(BUILD_DIR, "input_gen.txt")
    with open(filename, 'w') as f:
        f.writelines(output_lines)
    return filename

def main():
    if not compile_renderer(): return
    
    print(f"🚀 Starting Animation ({TOTAL_FRAMES} frames) - ORIGINAL VERSION")
    print(f"⚠️  Warning: Without BVH, this might be very slow.")
    
    for i in range(TOTAL_FRAMES):
        print(f"   [Frame {i+1}/{TOTAL_FRAMES}] Rendering...")
        input_file = generate_frame_input(i)
        if not input_file: break
        
        # 執行
        cmd = [os.path.abspath(EXE_PATH), "input_gen.txt"]
        
        start_t = time.time()
        # 顯示 stdout 讓你確認它還活著
        p = subprocess.run(cmd, cwd=BUILD_DIR) 
        end_t = time.time()
        
        if p.returncode != 0:
            print(f"❌ Error rendering frame {i}")
            break
            
        print(f"      -> Time: {end_t - start_t:.2f}s")

        # 移動圖片
        src_ppm = os.path.join(BUILD_DIR, "AdvCG_mis.ppm")
        target_ppm = os.path.join(OUTPUT_DIR, f"frame_{i:03d}.ppm")
        
        if os.path.exists(src_ppm):
            shutil.move(src_ppm, target_ppm)
        else:
            print(f"⚠️ Warning: Output image not found for frame {i}")

    print(f"\n✅ Animation Done! Check '{OUTPUT_DIR}' folder.")

if __name__ == "__main__":
    main()