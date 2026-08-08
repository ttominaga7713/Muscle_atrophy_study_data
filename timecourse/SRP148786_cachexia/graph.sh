#!/bin/bash

# エラー発生時に停止
set -e

# 1. ヘルプオプション (-h または --help) のみの処理
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    echo "========================================================================"
    echo " 📊 グラフ一括作成スクリプト (graph.sh) 統合マニュアル"
    echo "========================================================================"
    echo "【基本的な使い方】"
    echo "  $ bash $0 <list.csv> <output.csv>"
    echo ""
    echo "  1. <list.csv>   : graph.sh に渡す全ファイルリスト (タイムポイント設定用)"
    echo "  2. <output.csv> : グラフ作成の指示書 (遺伝子リスト, タイトル, [左上ラベル])"
    echo ""
    echo "------------------------------------------------------------------------"
    echo "【1. 指示書 (output.csv) の書き方】"
    echo "  カンマ(,) または タブ(\t) 区切りで作成します。1行目のヘッダーはスキップされます。"
    echo "  ※ 3列目に (A) などの文字を入れると、グラフの左上に独立して表示されます（空欄でもOK）。"
    echo "  ※ 画像拡張子(.pngなど)は書く必要はありません！"
    echo "========================================================================"
    exit 0
fi

# 2. 引数が足りない場合の短いエラー表示
if [ "$#" -ne 2 ]; then
    echo "⚠️ エラー: 引数が不足しています。"
    echo "使用方法: bash $0 <list.csv> <output.csv>"
    echo "💡 詳しいマニュアルを見るには bash $0 -h を実行してください。"
    exit 1
fi

LIST_CSV=$1
OUTPUT_CSV=$2

echo "📊 指示書 ($OUTPUT_CSV) に基づきグラフの一括作成を開始します..."

# ファイルディスクリプタ9を使用してループ
while IFS= read -r -u 9 line; do
    # 空行をスキップ
    if [[ -z "$line" ]]; then continue; fi
    
    # Windows特有の改行コード(\r)を削除
    line=$(echo "$line" | tr -d '\r')
    
    # ヘッダー行 (gene filename label) をスキップ
    if [[ "$line" =~ ^gene[[:space:],]*filename ]]; then continue; fi

    # タブ区切りかカンマ区切りかを自動判定して3列分を分割
    if [[ "$line" == *$'\t'* ]]; then
        gene_csv=$(echo "$line" | awk -F'\t' '{print $1}')
        filename_raw=$(echo "$line" | awk -F'\t' '{print $2}')
        label_raw=$(echo "$line" | awk -F'\t' '{print $3}')
    else
        gene_csv=$(echo "$line" | awk -F',' '{print $1}')
        filename_raw=$(echo "$line" | awk -F',' '{print $2}')
        label_raw=$(echo "$line" | awk -F',' '{print $3}')
    fi
    
    # 余計なダブルクォーテーションや前後の空白を削除
    gene_csv=$(echo "$gene_csv" | tr -d '"' | xargs)
    filename_raw=$(echo "$filename_raw" | tr -d '"' | xargs)
    label_text=$(echo "$label_raw" | tr -d '"' | xargs)

    # 拡張子を消しつつ、末尾に紛れ込んだカンマ(,)や空白も強制的に削除する
    title=$(echo "$filename_raw" | sed 's/\.[a-zA-Z0-9]*$//' | sed 's/[,[:space:]]*$//')

    echo "=================================================="
    if [ -n "$label_text" ]; then
        echo "📈 処理中: $title (パネルラベル: $label_text)"
    else
        echo "📈 処理中: $title"
    fi
    echo "📁 使用リスト: $gene_csv"
    
    # 指定された遺伝子リストのファイルが存在するかチェック
    if [ ! -f "$gene_csv" ]; then
        echo "⚠️ エラー: '$gene_csv' が見つかりません。スキップします。"
        continue
    fi

    # Python/Docker処理
    docker run -i --rm \
      -v "$(pwd):/work" \
      -w /work \
      ezojika7713/rnaseq-meta:latest \
      python3 - "$LIST_CSV" "$gene_csv" "$title" "$label_text" << 'EOF'

import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import os
import warnings
import re
import math

warnings.filterwarnings('ignore')

FILE_CSV = sys.argv[1]
GENE_CSV = sys.argv[2]
GRAPH_TITLE = sys.argv[3]
# 4番目の引数（ラベル）があれば受け取り、なければ空文字にする
LABEL_TEXT = sys.argv[4] if len(sys.argv) > 4 else ""

safe_filename = re.sub(r'[\\/:\*\?"<>\|]', '_', GRAPH_TITLE)
OUTPUT_PLOT = f"{safe_filename}.png"

def get_optimal_ticks(max_data, max_asterisk):
    if max_data <= 1.2:
        step = 0.2
        final_upper = 1.2
        if max_asterisk > final_upper * 0.95:
            final_upper += step
        return np.arange(0, final_upper + 0.01, step), final_upper

    if max_data <= 2.0:
        step = 0.2
    elif max_data <= 5:
        step = 1
    elif max_data <= 10:
        step = 2
    elif max_data <= 20:
        step = 5
    elif max_data <= 50:
        step = 10
    else:
        step = 20
        
    normal_upper = math.ceil(max_data / step) * step
    
    threshold = normal_upper - (step / 2.0)
    if max_data >= threshold:
        final_upper = normal_upper + step
    else:
        final_upper = normal_upper
        
    if max_asterisk > final_upper * 0.95:
        final_upper += step
        
    return np.arange(0, final_upper + (step * 0.1), step), final_upper

def main():
    genes_df = pd.read_csv(GENE_CSV)
    files_df = pd.read_csv(FILE_CSV)
    files_df['timepoint'] = files_df['timepoint'].astype(str)
    timepoints_order = files_df['timepoint'].tolist()
    
    all_dfs = {}
    constants = {}
    
    for index, row in files_df.iterrows():
        fname = str(row['filename'])
        tp = str(row['timepoint'])
        
        if os.path.exists(fname):
            df = pd.read_csv(fname, sep='\t')
            if 'logFC' in df.columns:
                df['FC'] = 2 ** df['logFC']
            all_dfs[tp] = df
        else:
            try:
                val = float(fname)
                constants[tp] = val
            except ValueError:
                print(f"Warning: '{fname}' はファイルとして存在せず、数値でもありません。スキップします。")

    fig = plt.figure(figsize=(11.5, 6))
    ax1 = fig.add_axes([0.10, 0.15, 0.59, 0.70])

    ax2 = None
    if genes_df.shape[1] >= 2:
        ax2 = ax1.twinx()

    x_positions = range(len(timepoints_order))
    lines = []
    labels = []

    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_idx = 0

    ax1_max_data = 0.0
    ax1_max_asterisk = 0.0

    for col_idx in range(min(genes_df.shape[1], 2)):
        target_ax = ax1 if col_idx == 0 else ax2
        linestyle = '-' if col_idx == 0 else '--'
        col_name = genes_df.columns[col_idx]
        
        for gene_id in genes_df[col_name].dropna():
            fc_values = []
            fdr_values = []
            symbol = gene_id
            
            for tp in timepoints_order:
                if tp in all_dfs:
                    row = all_dfs[tp][all_dfs[tp]['Geneid'] == gene_id]
                    if not row.empty:
                        fc_values.append(row.iloc[0]['FC'])
                        fdr_values.append(row.iloc[0]['FDR'])
                        symbol = row.iloc[0].get('GeneSymbol', gene_id)
                    else:
                        fc_values.append(np.nan)
                        fdr_values.append(np.nan)
                elif tp in constants:
                    fc_values.append(constants[tp])
                    fdr_values.append(np.nan)
                else:
                    fc_values.append(np.nan)
                    fdr_values.append(np.nan)

            if col_idx == 0:
                for fc, fdr in zip(fc_values, fdr_values):
                    if pd.notna(fc):
                        ax1_max_data = max(ax1_max_data, fc)
                        if pd.notna(fdr) and fdr < 0.05:
                            ax1_max_asterisk = max(ax1_max_asterisk, fc * 1.15)

            current_color = color_cycle[color_idx % len(color_cycle)]
            color_idx += 1

            ln, = target_ax.plot(x_positions, fc_values, marker='o', linestyle=linestyle, color=current_color, label=symbol, linewidth=2)
            lines.append(ln)
            labels.append(symbol)

            for i, (fc, fdr) in enumerate(zip(fc_values, fdr_values)):
                if pd.notna(fdr) and pd.notna(fc):
                    mark = "**" if fdr < 0.01 else "*" if fdr < 0.05 else ""
                    if mark:
                        target_ax.text(i, fc + (fc * 0.08), mark, ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax1.set_xlabel('Timepoint (Days)', fontsize=22)
    ax1.set_ylabel('Fold Change', fontsize=22)
    
    ax1.yaxis.set_label_coords(-0.09, 0.5)
    
    ax1_ticks, ax1_upper = get_optimal_ticks(ax1_max_data, ax1_max_asterisk)
    ax1.set_yticks(ax1_ticks)
    ax1.set_ylim(0, ax1_upper)
    
    if ax2:
        ax2.margins(y=0.2)
    
    ax1.tick_params(axis='both', which='major', labelsize=18, direction='out')
    if ax2:
        ax2.tick_params(axis='y', which='major', labelsize=18, direction='out')
    
    ax1.set_xticks(x_positions)
    ax1.set_xticklabels(timepoints_order)
    
    ax1.axhline(1, color='gray', linestyle=':', linewidth=1)
    
    ax1.legend(lines, labels, bbox_to_anchor=(1.08, 1.0), loc='upper left', fontsize=16, handlelength=3)
    
    # 【変更箇所】タイトルのフォントサイズを28から24へ若干小さくしました
    fig.suptitle(GRAPH_TITLE, fontsize=24, y=0.95)
    
    if LABEL_TEXT:
        fig.text(0.02, 0.98, LABEL_TEXT, fontsize=28, va='top', ha='left')
    
    ax1.set_ylim(bottom=0)
    
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"✅ グラフを {OUTPUT_PLOT} に保存しました。")

if __name__ == "__main__":
    main()
EOF

done 9< "$OUTPUT_CSV"

echo "=================================================="
echo "🎉 すべてのグラフの一括作成が完了しました！"