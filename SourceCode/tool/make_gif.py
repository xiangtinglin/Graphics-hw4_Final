'''
輸出video:
ffmpeg -framerate 24 -i animation_frames_original/frame_%03d.ppm -c:v libx264 -pix_fmt yuv420p animation_frames_original/output.mp4
'''
import os
import glob
from PIL import Image

# ================= 設定 =================
INPUT_DIR = "animation_frames_original"        # 圖片所在的資料夾
OUTPUT_GIF = "result_animation.gif"   # 輸出的 GIF 檔名
OUTPUT_DIR = "animation_frames_original"
FPS = 24                              # 每秒幾張 (24 為電影標準)
DURATION = int(1000 / FPS)            # 每張圖停留毫秒數

def create_gif():
    # 1. 搜尋並排序圖片
    # 確保路徑正確，例如 animation_frames/frame_000.ppm
    search_pattern = os.path.join(INPUT_DIR, "frame_*.ppm")
    files = sorted(glob.glob(search_pattern))
    
    if not files:
        print(f"❌ 錯誤：在 '{INPUT_DIR}' 找不到任何 ppm 圖片。")
        print("   請確認您是否已經執行過 render_animation 腳本，且圖片確實存在。")
        return

    print(f"📂 找到 {len(files)} 張圖片，正在合成 GIF...")

    images = []
    for filename in files:
        try:
            # 開啟圖片
            img = Image.open(filename)
            # PPM 有時需要轉為 RGB 模式才能正確儲存為 GIF
            img = img.convert('RGB')
            # 為了讓 GIF 檔案小一點，可以考慮縮圖 (選擇性)
            # img.thumbnail((512, 512)) 
            images.append(img)
            print(f"   讀取: {filename}", end='\r')
        except Exception as e:
            print(f"⚠️ 無法讀取 {filename}: {e}")

    print("\n💾 正在儲存 GIF (這可能需要幾秒鐘)...")
    
    if images:
        # 2. 儲存為 GIF
        images[0].save(
            os.path.join(OUTPUT_DIR,OUTPUT_GIF),
            save_all=True,
            append_images=images[1:],
            optimize=True,   # 開啟最佳化以縮小檔案
            duration=DURATION,
            loop=0           # 0 代表無限循環
        )
        print(f"🎉 成功！動畫已儲存為：{os.path.join(OUTPUT_DIR,OUTPUT_GIF)}")
    else:
        print("❌ 沒有有效的圖片可以合成。")

if __name__ == "__main__":
    # 檢查套件
    try:
        import PIL
        create_gif()
    except ImportError:
        print("❌ 缺少 'Pillow' 套件。")
        print("請執行: pip install Pillow")