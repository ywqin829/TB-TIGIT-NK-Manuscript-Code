# ======================================================================
# AG分析SE.py — AlphaGenome-based Super-Enhancer Analysis
# Predicts H3K27ac signal and calls super-enhancers via AlphaGenome API
# ======================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from alphagenome.models import dna_client
from alphagenome.data import genome
import os

# ======================================================================
# Configuration
# ======================================================================

# TODO: Set your AlphaGenome API key here
MY_API_KEY = "YOUR_API_KEY_HERE"

CELL_TYPE = "CL:0000623"

# Output directory (relative to working directory)
SAVE_PATH = os.path.join(os.getcwd(), "output")
os.makedirs(SAVE_PATH, exist_ok=True)

gene_info = {
    "name": "TBX21",
    "chr": "chr17",
    "tss": 47733236,
    "strand": "+",
    "window_size": 524288  # Total width = 1,048,576 bp
}

print("===== AlphaGenome Super-Enhancer Analysis =====")
print(f"Output path: {SAVE_PATH}")
print("Starting analysis...")

# ======================================================================
# Core analysis function
# ======================================================================

def analyze_enhancer_intensity(gene_dict, api_key, cell_type):
    """Predict H3K27ac signal and call enhancer peaks via AlphaGenome."""

    print(f"\n--- Analyzing gene: {gene_dict['name']} ---")

    client = dna_client.create(api_key)

    interval = genome.Interval(
        chromosome=gene_dict["chr"],
        start=gene_dict["tss"] - gene_dict["window_size"],
        end=gene_dict["tss"] + gene_dict["window_size"]
    )

    print(f"Requesting region: {interval.chromosome}:{interval.start}-{interval.end}")
    print(f"Window width: {interval.end - interval.start} bp")

    # Request H3K27ac prediction from AlphaGenome
    output = client.predict_interval(
        interval=interval,
        requested_outputs=[dna_client.OutputType.CHIP_HISTONE],
        ontology_terms=[cell_type],
    )

    # Extract H3K27ac track
    h3k27ac_meta = output.chip_histone.metadata[
        output.chip_histone.metadata['histone_mark'] == 'H3K27ac'
    ]

    if h3k27ac_meta.empty:
        print("ERROR: No H3K27ac data found")
        return None, None, None

    idx = h3k27ac_meta.index[0]
    signal_low = output.chip_histone.values[:, idx]

    # Upsample to base-pair resolution
    positions = np.arange(interval.start, interval.end)
    res_factor = len(positions) // len(signal_low)
    signal = np.repeat(signal_low, res_factor)
    print(f"Data retrieved successfully. Signal length: {len(signal)}")

    # Peak detection (mean + 2 SD threshold)
    threshold = np.mean(signal) + 2 * np.std(signal)
    above = signal > threshold
    diff = np.diff(above.astype(int))

    starts = np.where(diff == 1)[0] + 1
    ends = np.where(diff == -1)[0] + 1

    if above[0]:
        starts = np.insert(starts, 0, 0)
    if above[-1]:
        ends = np.append(ends, len(above) - 1)

    # Score each peak
    results = []
    for i, (s, e) in enumerate(zip(starts, ends)):
        s_pos = interval.start + s
        e_pos = interval.start + e
        sub_sig = signal[s:e]
        max_val = np.max(sub_sig)
        width = e_pos - s_pos
        score = max_val * np.log10(width)

        results.append({
            "ID": i + 1,
            "Chromosome": interval.chromosome,
            "Start": s_pos,
            "End": e_pos,
            "Width": width,
            "Max_Val": max_val,
            "SE_Score": score
        })

    df = pd.DataFrame(results)
    df['Is_SE'] = df['SE_Score'] >= np.percentile(df['SE_Score'], 75)

    return df, signal, positions, interval

# ======================================================================
# Run analysis
# ======================================================================

print("\n===== Running enhancer analysis =====")

df_intensity, h3k27ac_signal, positions, interval = analyze_enhancer_intensity(
    gene_info, MY_API_KEY, CELL_TYPE
)

# ======================================================================
# Save CSV results
# ======================================================================

if df_intensity is not None:

    print("\n===== Saving CSV results =====")

    # Filter super-enhancers
    super_enhancers = df_intensity[df_intensity['Is_SE']].copy()
    super_enhancers = super_enhancers.sort_values(
        'SE_Score', ascending=False
    ).reset_index(drop=True)

    # Add annotation columns
    super_enhancers['Width_kb'] = super_enhancers['Width'] / 1000
    super_enhancers['Distance_to_TSS'] = super_enhancers['Start'] - gene_info['tss']
    super_enhancers['Direction'] = super_enhancers['Distance_to_TSS'].apply(
        lambda x: 'Upstream' if x < 0 else 'Downstream'
    )

    # All enhancers
    csv_all_path = os.path.join(SAVE_PATH, "TBX21_All_Enhancers.csv")
    df_intensity.to_csv(csv_all_path, index=False)
    print(f"All enhancers saved: {csv_all_path}")

    # Super-enhancers only
    csv_se_path = os.path.join(SAVE_PATH, "TBX21_Super_Enhancers.csv")
    super_enhancers.to_csv(csv_se_path, index=False)
    print(f"Super-enhancers saved: {csv_se_path}")

    # Preview
    print("\nSuper-enhancer preview:")
    print("=" * 80)
    print(super_enhancers[['ID', 'Chromosome', 'Start', 'End', 'Width_kb',
                           'Max_Val', 'SE_Score', 'Distance_to_TSS']].to_string(index=False))

    print(f"\nTotal enhancers detected: {len(df_intensity)}")
    print(f"Super-enhancers: {len(super_enhancers)}")

# ======================================================================
# Visualization (PDF)
# ======================================================================

    print("\n===== Generating visualization (PDF) =====")

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

    # --- Panel A: SE ranking bar chart ---
    ax1 = fig.add_subplot(gs[0, :])

    df_se_sorted = super_enhancers.sort_values(
        'SE_Score', ascending=False
    ).reset_index(drop=True)
    df_se_sorted['Rank'] = range(1, len(df_se_sorted) + 1)

    colors = ['#FFD700' if i < 3 else '#FFA500'
              for i in range(len(df_se_sorted))]
    bars = ax1.bar(df_se_sorted['ID'], df_se_sorted['SE_Score'],
                   color=colors, edgecolor='black', linewidth=1.5)

    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2., height + 100,
                 f"{int(height)}", ha='center', va='bottom', fontsize=9)

    ax1.set_xlabel('Super-Enhancer ID', fontsize=12, fontweight='bold')
    ax1.set_ylabel('SE Score', fontsize=12, fontweight='bold')
    ax1.set_title('Panel A: TBX21 Super-Enhancers Ranking by SE Score',
                  fontsize=13, fontweight='bold')
    ax1.grid(True, axis='y', alpha=0.3)

    # --- Panel B: width vs intensity scatter ---
    ax2 = fig.add_subplot(gs[1, 0])

    se_data = df_intensity[df_intensity['Is_SE']]
    typical_data = df_intensity[~df_intensity['Is_SE']]

    ax2.scatter(typical_data['Width'] / 1000, typical_data['Max_Val'],
                c='skyblue', s=150, alpha=0.7, edgecolors='black', linewidth=1.5,
                label='Typical Enhancer')
    ax2.scatter(se_data['Width'] / 1000, se_data['Max_Val'],
                c='gold', s=200, alpha=0.8, edgecolors='black', linewidth=2,
                label='Super-Enhancer')

    for idx, row in se_data.head(3).iterrows():
        ax2.text(row['Width'] / 1000, row['Max_Val'] * 1.05,
                 f"SE_{row['ID']}", ha='center', fontsize=9, fontweight='bold')

    ax2.set_xlabel('Enhancer Width (kb)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Max Signal Intensity', fontsize=11, fontweight='bold')
    ax2.set_title('Panel B: SE Characteristics', fontsize=12, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3)

    # --- Panel C: distance to TSS distribution ---
    ax3 = fig.add_subplot(gs[1, 1])

    df_sorted_dist = super_enhancers.sort_values('Distance_to_TSS')
    colors_dist = ['green' if d < 0 else 'red'
                   for d in df_sorted_dist['Distance_to_TSS']]

    bars_dist = ax3.barh(range(len(df_sorted_dist)),
                         df_sorted_dist['Distance_to_TSS'] / 1000,
                         color=colors_dist, alpha=0.7,
                         edgecolor='black', linewidth=1.5)

    for i, (idx, row) in enumerate(df_sorted_dist.iterrows()):
        direction = chr(8593) if row['Distance_to_TSS'] < 0 else chr(8595)
        ax3.text(row['Distance_to_TSS'] / 1000, i,
                 f"  {direction} {abs(row['Distance_to_TSS']) / 1000:.1f}kb",
                 va='center', fontsize=9, fontweight='bold')

    ax3.axvline(0, color='black', linestyle='--', linewidth=2, alpha=0.5)
    ax3.set_yticks(range(len(df_sorted_dist)))
    ax3.set_yticklabels([f"SE_{row['ID']}"
                         for _, row in df_sorted_dist.iterrows()])
    ax3.set_xlabel('Distance from TSS (kb)', fontsize=11, fontweight='bold')
    ax3.set_title('Panel C: SE Distribution Relative to TSS',
                  fontsize=12, fontweight='bold')
    ax3.grid(True, axis='x', alpha=0.3)

    # --- Panel D: genomic signal track ---
    ax4 = fig.add_subplot(gs[2, :])

    ax4.plot(positions, h3k27ac_signal, color='#2E86AB', linewidth=1.5,
             label='H3K27ac (Predicted)', alpha=0.8)
    ax4.fill_between(positions, h3k27ac_signal, color='#2E86AB', alpha=0.3)

    ax4.axvline(gene_info['tss'], color='red', linestyle='--', linewidth=3,
                label='TSS', zorder=10)
    ax4.text(gene_info['tss'], ax4.get_ylim()[1] * 0.92,
             '  TBX21 TSS', color='red', fontweight='bold',
             fontsize=11, ha='center')

    for idx, row in super_enhancers.iterrows():
        ax4.axvspan(row['Start'], row['End'], color='gold', alpha=0.4,
                    label='Super-Enhancer' if idx == 0 else "")
        mid_pos = (row['Start'] + row['End']) // 2
        ax4.text(mid_pos, ax4.get_ylim()[1] * 0.85,
                 f"SE_{row['ID']}", ha='center', fontsize=9,
                 fontweight='bold', color='#8B4513')

    ax4.set_xlabel('Genomic Position (chr17)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('H3K27ac Signal Intensity', fontsize=12, fontweight='bold')
    ax4.set_title('Panel D: TBX21 Genomic Landscape with Super-Enhancers',
                  fontsize=13, fontweight='bold')
    ax4.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax4.grid(True, alpha=0.3)

    fig.suptitle(
        f'TBX21 Super-Enhancers Analysis\n'
        f'(Total: {len(df_intensity)} enhancers, '
        f'{len(super_enhancers)} super-enhancers)',
        fontsize=15, fontweight='bold', y=0.99
    )

    pdf_path = os.path.join(SAVE_PATH, "TBX21_Super_Enhancers_Visualization.pdf")
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight', dpi=300)
    print(f"\nPDF saved: {pdf_path}")
    plt.close()

# ======================================================================
# Summary report (TXT)
# ======================================================================

    print("\n===== Generating summary report =====")

    summary = f"""
{'=' * 80}
TBX21 Super-Enhancer Analysis Summary
{'=' * 80}

Analysis region: {interval.chromosome}:{interval.start:,}-{interval.end:,}
TSS position: {gene_info['chr']}:{gene_info['tss']:,}

Total enhancers detected: {len(df_intensity)}
Super-enhancers: {len(super_enhancers)}
Super-enhancer ratio: {len(super_enhancers) / len(df_intensity) * 100:.1f}%

Top 3 super-enhancers:
"""
    for i, (idx, row) in enumerate(super_enhancers.head(3).iterrows(), 1):
        summary += f"""
  {i}. SE_{row['ID']}
     - Location: {row['Chromosome']}:{row['Start']:,}-{row['End']:,}
     - Width: {row['Width_kb']:.1f} kb
     - SE Score: {row['SE_Score']:.2f}
     - Distance to TSS: {abs(row['Distance_to_TSS']) / 1000:.1f} kb ({row['Direction']})
"""

    summary += f"\n{'=' * 80}\n"
    print(summary)

    txt_path = os.path.join(SAVE_PATH, "TBX21_Analysis_Summary.txt")
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write(summary)
    print(f"Summary saved: {txt_path}")

else:
    print("ERROR: Analysis failed")

print("\n===== AG分析SE.py complete =====")
print(f"All outputs in: {SAVE_PATH}")
