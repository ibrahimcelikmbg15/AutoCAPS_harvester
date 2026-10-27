# -*- coding: utf-8 -*-
import streamlit as st
import pandas as pd
import re
import io
import matplotlib.subplots as plt_sub
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import plotly.graph_objects as go

st.set_page_config(page_title="AutoCAPS Harvester", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
    <style>
    .metric-box {
        background-color: rgba(14, 165, 233, 0.1);
        border: 1px solid #0ea5e9;
        border-radius: 8px;
        padding: 15px;
        text-align: center;
        margin-bottom: 10px;
    }
    .metric-value { font-size: 2rem; font-weight: bold; color: #0ea5e9; }
    .metric-label { font-size: 1rem; font-weight: 600; }
    .credit-text {
        font-size: 1.1rem;
        color: #10b981;
        font-weight: 700;
        border-left: 4px solid #10b981;
        padding-left: 15px;
        margin-top: 5px;
    }
    </style>
""", unsafe_allow_html=True)

st.markdown("<h1>🌾 AutoCAPS Harvester: Interactive Results Viewer</h1>", unsafe_allow_html=True)
st.markdown("### Cloud-Based Dynamic Dashboard for Pre-Processed CAPS Markers")
st.markdown("<div class='credit-text'>Developed by: Assoc. Prof. Dr. Ibrahim Celik</div>", unsafe_allow_html=True)
st.divider()

def plot_vertical_map(df, chrom_lengths, selected_enzyme, sel_chroms):
    fig = go.Figure()
    
    if selected_enzyme == "All Successful Enzymes":
        plot_df = df[df['Compatible_Enzyme_Count'] > 0]
    else:
        plot_df = df[df['Used_Enzyme_List'].apply(lambda x: isinstance(x, list) and selected_enzyme in x)]

    # Sadece seçilen kromozomları listele (boşsa hepsini al)
    if not sel_chroms:
        valid_chroms = sorted([ch for ch in chrom_lengths.keys() if chrom_lengths[ch] > 0])
    else:
        valid_chroms = sorted([ch for ch in sel_chroms if chrom_lengths.get(ch, 0) > 0])
    
    max_len_mb = 0
    for chrom in valid_chroms:
        length_mb = chrom_lengths[chrom] / 1_000_000
        if length_mb > max_len_mb:
            max_len_mb = length_mb
        
        # Kromozom çubuğunu çiz (Sadece seçilenler çizilir)
        fig.add_trace(go.Scatter(
            x=[chrom, chrom], y=[0, length_mb], 
            mode='lines', line=dict(color='lightgray', width=8), 
            showlegend=False, hoverinfo='skip'
        ))
        
        chrom_df = plot_df[plot_df['FASTA_Header'] == chrom]
        if not chrom_df.empty:
            y_vals = chrom_df['POS'] / 1_000_000
            hover_texts = chrom_df.apply(lambda r: f"SNP: {r['SNP_ID']}<br>POS: {r['POS']}<br>Enzymes: {len(r.get('Used_Enzyme_List', []))}", axis=1)
            
            fig.add_trace(go.Scatter(
                x=[chrom] * len(y_vals), y=y_vals,
                mode='markers',
                marker=dict(symbol='line-ew', color='#0ea5e9', size=25, line=dict(width=3, color='#0ea5e9')),
                text=hover_texts, hoverinfo='text', showlegend=False
            ))

    fig.update_layout(
        title=f"Interactive Physical Map: {selected_enzyme}",
        yaxis_title="Physical Position (Mb)",
        yaxis=dict(autorange='reversed', range=[max_len_mb * 1.05, - (max_len_mb * 0.05)]),
        xaxis_title="Chromosomes",
        plot_bgcolor='rgba(0,0,0,0)',
        height=700,
        margin=dict(l=40, r=40, t=60, b=40)
    )
    return fig

def parse_profile_for_gel(profile_str, enzyme, ref, alt, pcr_size):
    frag_data = {ref: [pcr_size], alt: [pcr_size]}
    if pd.isna(profile_str) or not str(profile_str).strip() or "No Suitable" in str(profile_str):
        return frag_data
    
    enzymes_profiles = str(profile_str).split(" || ")
    target_prof = None
    for prof in enzymes_profiles:
        if prof.startswith(enzyme + " ->"):
            target_prof = prof
            break
            
    if not target_prof: return frag_data
    
    alleles_split = target_prof.split(" | ")
    for part in alleles_split:
        m = re.search(r'Allele (.*?): \[(.*?)\]', part)
        if m:
            al = m.group(1).strip()
            sizes_str = m.group(2)
            if "No Cut" in sizes_str:
                frag_data[al] = [pcr_size]
            else:
                sizes = [int(s.replace("bp", "").strip()) for s in sizes_str.split("+")]
                frag_data[al] = sizes
    return frag_data

def render_gel_simulator(frag_data, ref, alt, pcr_size, enzyme):
    fig, ax = plt.subplots(figsize=(6, 8))
    
    max_bp = ((pcr_size // 100) + 1) * 100 if pcr_size % 100 != 0 else pcr_size
    marker_bands = sorted(list(range(50, max_bp + 100, 50)), reverse=True)
    
    frags1 = frag_data.get(ref, [pcr_size])
    frags2 = frag_data.get(alt, [pcr_size])
    
    lanes = {
        "DNA\nMarker": marker_bands,
        f"Allele 1\n({ref}{ref})": frags1,
        f"Allele 2\n({alt}{alt})": frags2,
        f"Hybrid\n({ref}{alt})": sorted(list(set(frags1 + frags2)), reverse=True)
    }
    
    well_y = max_bp + 50
    for i in range(1, len(lanes) + 1):
        well = patches.Rectangle((i - 0.25, well_y), 0.5, max_bp * 0.03, fill=False, edgecolor='black', lw=2)
        ax.add_patch(well)
        ax.text(i, well_y + (max_bp * 0.04), "WELL", color='gray', ha='center', fontsize=9, fontweight='bold')
        
    for i, (lane_name, bands) in enumerate(lanes.items()):
        x = i + 1
        for band in bands:
            alpha_core = 0.9 if lane_name != "DNA\nMarker" else 0.5
            b_color = '#0ea5e9' if lane_name != "DNA\nMarker" else 'gray'
            ax.hlines(band, x - 0.25, x + 0.25, colors=b_color, linewidth=18, alpha=alpha_core * 0.3)
            ax.hlines(band, x - 0.2, x + 0.2, colors='black', linewidth=4, alpha=alpha_core) 
            if lane_name == "DNA\nMarker":
                ax.text(x - 0.35, band, f"{band}", color='black', va='center', ha='right', fontsize=10)

    ax.set_xticks(list(range(1, len(lanes) + 1)))
    ax.set_xticklabels(lanes.keys(), fontsize=11, fontweight='bold')
    ax.set_ylim(0, max_bp + (max_bp * 0.15))
    ax.set_ylabel("DNA Size (bp)", fontsize=12, fontweight='bold')
    ax.set_title(f"{enzyme} Restriction Profile", color='#0ea5e9', fontweight='bold', fontsize=14, pad=20)
    return fig

# --- SIDEBAR: DATA UPLOAD ---
st.sidebar.header("📂 Upload Results")
uploaded_file = st.sidebar.file_uploader("Upload AutoCAPS Report (.xlsx or .csv)", type=["xlsx", "csv"])

if uploaded_file:
    try:
        if uploaded_file.name.endswith('.csv'):
            res_df = pd.read_csv(uploaded_file)
        else:
            res_df = pd.read_excel(uploaded_file)
            
        if 'Used_Enzyme_List' not in res_df.columns:
            def extract_enzymes_from_profile(profile_text):
                if pd.isna(profile_text) or "No Suitable" in str(profile_text) or "ERROR" in str(profile_text):
                    return []
                return [p.split(" ->")[0].strip() for p in str(profile_text).split(" || ") if " ->" in p]
            
            res_df['Used_Enzyme_List'] = res_df['CAPS_Restriction_Profile'].apply(extract_enzymes_from_profile)
        else:
            res_df['Used_Enzyme_List'] = res_df['Used_Enzyme_List'].apply(
                lambda x: [e.strip() for e in str(x).split(',')] if pd.notna(x) and str(x).strip() else []
            )
        
        all_enzymes = []
        for enz_list in res_df['Used_Enzyme_List']:
            if isinstance(enz_list, list):
                all_enzymes.extend(enz_list)
        unique_enzymes = sorted(list(set(all_enzymes)))
        
        chrom_lengths = {}
        for ch in res_df['FASTA_Header'].unique():
            if ch != "-" and pd.notna(ch):
                max_pos = res_df[res_df['FASTA_Header'] == ch]['POS'].max()
                chrom_lengths[ch] = int(max_pos * 1.05) if pd.notna(max_pos) else 1000000

        # --- DASHBOARD UI ---
        st.markdown("## 🔍 Dynamic Filtering")
        c1, c2, c3, c4 = st.columns(4)
        
        # Çoklu seçim kutusuna geçiş
        all_chroms = sorted([ch for ch in res_df['FASTA_Header'].unique() if ch != "-" and pd.notna(ch)])
        sel_chroms = c1.multiselect("Target Chromosome(s)", all_chroms, default=[], help="Leave empty to select all chromosomes")
        
        sel_start = c2.number_input("Region Start (bp)", value=0, min_value=0)
        sel_end = c3.number_input("Region End (bp) (0 = All)", value=0, min_value=0)
        
        all_enz_opts = ["All Successful Enzymes"] + unique_enzymes
        sel_enz = c4.selectbox("Restriction Enzyme", all_enz_opts)
        
        filt_df = res_df.copy()
        
        # Seçili kromozom varsa ona göre filtrele
        if sel_chroms:
            filt_df = filt_df[filt_df['FASTA_Header'].isin(sel_chroms)]
            
        if sel_end > 0:
            filt_df = filt_df[(filt_df['POS'] >= sel_start) & (filt_df['POS'] <= sel_end)]
        
        if sel_enz != "All Successful Enzymes":
            filt_df = filt_df[filt_df['Used_Enzyme_List'].apply(lambda x: isinstance(x, list) and sel_enz in x)]
            
        st.divider()

        st.markdown("## 📊 1. Marker Statistics")
        
        total_snp = len(filt_df)
        caps_success = len(filt_df[filt_df['Compatible_Enzyme_Count'] > 0])
        
        m1, m2 = st.columns(2)
        m1.markdown(f"<div class='metric-box'><div class='metric-value'>{total_snp}</div><div class='metric-label'>Filtered Markers in View</div></div>", unsafe_allow_html=True)
        m2.markdown(f"<div class='metric-box'><div class='metric-value'>{caps_success}</div><div class='metric-label'>Successful CAPS Conversions</div></div>", unsafe_allow_html=True)
        
        st.divider()
        
        st.markdown("## 🌍 2. Interactive Physical Map")
        if caps_success > 0:
            # Map fonksiyonuna seçilen kromozomları gönderiyoruz
            fig_map = plot_vertical_map(filt_df, chrom_lengths, sel_enz, sel_chroms)
            st.plotly_chart(fig_map, use_container_width=True)
        else:
            st.info("No successful CAPS markers found to map with current filters.")
            
        st.divider()
        
        st.markdown("## 🖥️ 3. Virtual Gel Simulator")
        gel_candidates = filt_df[filt_df['Compatible_Enzyme_Count'] > 0]
        
        if gel_candidates.empty:
            st.warning("No suitable polymorphic cut profile found for gel simulation.")
        else:
            col_gel1, col_gel2 = st.columns([1, 2])
            with col_gel1:
                selected_snp_id = st.selectbox("Select Target SNP:", gel_candidates['SNP_ID'].tolist())
                snp_data = gel_candidates[gel_candidates['SNP_ID'] == selected_snp_id].iloc[0]
                
                ref_allele = snp_data['REF']
                alt_allele = snp_data['ALT']
                pcr_size = snp_data['PCR_Size_bp']
                enz_list = snp_data.get('Used_Enzyme_List', [])
                profile_str = snp_data['CAPS_Restriction_Profile']
                
                selected_snp_enz = st.selectbox("Select Restriction Enzyme:", enz_list if isinstance(enz_list, list) and len(enz_list) > 0 else ["No Enzyme"])
                frag_dict = parse_profile_for_gel(profile_str, selected_snp_enz, ref_allele, alt_allele, pcr_size)
                
                st.info(f"**Target SNP:** {selected_snp_id}\n\n"
                        f"**PCR Amplicon Size:** {pcr_size} bp\n\n"
                        f"**Allele 1 ({ref_allele}):** {frag_dict.get(ref_allele, [pcr_size])} bp\n\n"
                        f"**Allele 2 ({alt_allele}):** {frag_dict.get(alt_allele, [pcr_size])} bp")
                        
            with col_gel2:
                if selected_snp_enz != "No Enzyme":
                    fig_gel = render_gel_simulator(frag_dict, ref_allele, alt_allele, pcr_size, selected_snp_enz)
                    st.pyplot(fig_gel)
                
        st.divider()
        st.markdown("## 💾 4. Data Preview & Excel Export")
        st.dataframe(filt_df.drop(columns=['Used_Enzyme_List'], errors='ignore'))
        
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            export_df = filt_df.copy()
            if 'Used_Enzyme_List' in export_df.columns:
                export_df['Used_Enzyme_List'] = export_df['Used_Enzyme_List'].apply(
                    lambda x: ", ".join(x) if isinstance(x, list) else str(x)
                )
            export_df.to_excel(writer, sheet_name='Filtered_Markers', index=False)
            
        # Dosya ismini akıllı belirle
        chrom_name_tag = "Selected" if sel_chroms else "All"
        st.download_button(
            label="📂 Download Filtered Results as Excel (.xlsx)",
            data=output.getvalue(),
            file_name=f"AutoCAPS_Harvester_{chrom_name_tag}_Chromosomes_{sel_enz}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

    except Exception as e:
        st.error(f"Dosya işlenirken bir hata oluştu: {str(e)}")
else:
    st.info("👈 Please upload your `AutoCAPS_...csv` or `.xlsx` report from the sidebar to visualize the dashboard.")
