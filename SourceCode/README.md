# 目錄:
   ## Final output & 執行檔:
     * Final output資料夾:SourceCode/FinalOutput_and_exe
     * 執行檔:SourceCode/FinalOutput_and_exe/Final
     * v.1舊版output:SourceCode/Output_v.1
     
     * 測試編譯指令(先到根目錄./Main):
        clang++ -std=c++11 AdvCG_Final_MIS.cc -o ./output1206/Final
        ./output1206/Final input.txt
        
   ## BVH結果（Dashboard）:
     * 資料夾路徑:SourceCode/BVH
     * 匯出檔案:（1）BVH數據分析:SourceCode/BVH/benchmark_output/benchmark_report.html
               （2）數據資料:SourceCode/BVH/benchmark_output/benchmark_results.json
               （3）v.2的進階動圖都在:SourceCode/BVH/benchmark_output
     * 測試執行:（1）先到根目錄:SourceCode/BVH
               （2）再執行:python benchmark_bvh.py
