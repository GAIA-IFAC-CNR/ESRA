import streamlit as st
import pandas as pd
import numpy as np
import peptides
import io

st.set_page_config(page_title="ESRA - Epitope Selection Rational Approach", layout="wide")

st.title("ESRA - Epitope Selection Rational Approach")

# ==============================
# INPUT
# ==============================

name = st.text_input("Sequence name")
sequence = st.text_input("Sequence")

#st.subheader("Consensus Rules Thresholds")

rules_config = {
    'E4':               {'op': '<=', 'thresh': -0.091},
    'ProtFP4':          {'op': '>',  'thresh': 0.880},
    'SVGER4':           {'op': '<=', 'thresh': 0.195},
    'SVGER5':           {'op': '>',  'thresh': 0.416},
    'T3':               {'op': '<=', 'thresh': -0.113},
    'VHSE5':            {'op': '>',  'thresh': 0.011},
    'Z3':               {'op': '<=', 'thresh': -0.266},
    'charge':           {'op': '>',  'thresh': 1.542},
    'isoelectric_point':{'op': '>',  'thresh': 10.420},
    'PRIN3':            {'op': '>',  'thresh': 1.360}
}

# ==============================
# MAIN EXECUTION
# ==============================

col1, col2 = st.columns([8,1])
with col2:
    if st.button("Reset", type="primary"):

        # Svuota i campi input
        st.session_state["sequence"] = ""
        st.session_state["sequence_name"] = ""

        # Cancella risultati calcolati
        keys_to_clear = [
            "seq_df",
            "report_df",
            "styled_report",
            "name",
            "sequence"

        ]

        for key in keys_to_clear:
            if key in st.session_state:
                del st.session_state[key]

        st.rerun()



if st.button("Generate Report"):

    if not sequence:
        st.warning("Insert a sequence.")
        st.stop()

    # ==============================
    # 1. SUBSEQUENCE GENERATION
    # ==============================

    seq, first_index, last_index, leng, ids, c = [], [], [], [], [], []

    # Forward
    for start in range(len(sequence)):
        for i in range(7, 9):
            if (start + i) <= len(sequence):
                peptide = sequence[start:start+i]
                seq.append(peptide)
                first_index.append(start + 1)
                last_index.append(start + i)
                leng.append(i)
                ids.append(f'P{start+1}_L{i}')
                c.append('Y' if 'C' in peptide else 'N')

    # Reverse
    for start in range(len(sequence)-1, -1, -1):
        for i in range(7, 9):
            if (start - i) >= -1:
                peptide = sequence[start+1-i:start+1]
                if peptide not in seq:
                    seq.append(peptide)
                    first_index.append(start+1-i + 1)
                    last_index.append(start + 1)
                    leng.append(i)
                    ids.append(f'P{start+1-i + 1}_L{i}')
                    c.append('Y' if 'C' in peptide else 'N')

    # ==============================
    # 2. DESCRIPTOR CALCULATION
    # ==============================

    charge = []
    isoel = []
    E4 = []
    T3 = []
    VHSE5 = []
    SVGER5 = []
    ProtFP4 = []
    SVGER4 = []
    Z3 = []
    PRIN3 = []

    for pep in seq:
        peptide = peptides.Peptide(pep)

        charge.append(peptide.charge(pKscale="Bjellqvist"))
        isoel.append(peptide.isoelectric_point(pKscale="EMBOSS"))

        pcp = peptide.pcp_descriptors()
        E4.append(pcp[3])

        t_vals = peptide.t_scales()
        T3.append(t_vals[2])

        protfp = peptide.protfp_descriptors()
        ProtFP4.append(protfp[3])

        vhse_vals = peptide.vhse_scales()
        VHSE5.append(vhse_vals[4])

        svger_vals = peptide.svger_descriptors()
        SVGER5.append(svger_vals[4])
        SVGER4.append(svger_vals[3])

        z_vals = peptide.z_scales()
        Z3.append(z_vals[2])

        prin_vals = peptide.prin_components()
        PRIN3.append(prin_vals[2])

    seq_df = pd.DataFrame({
        'Sequence': seq,
        'First AA position': first_index,
        'Last AA position': last_index,
        'Length': leng,
        'ID': ids,
        'Cysteine': c,
        'charge': charge,
        'isoelectric_point': isoel,
        'E4': E4,
        'T3': T3,
        'VHSE5': VHSE5,
        'SVGER5': SVGER5,
        'SVGER4': SVGER4,
        'ProtFP4': ProtFP4,
        'Z3': Z3,
        'PRIN3': PRIN3
    })

    st.subheader("Subsequences and Descriptors")
    st.dataframe(seq_df, use_container_width=True)

    # ==============================
    # 3. CONSENSUS SCORE
    # ==============================

    def calculate_consensus(row):
        score = 0
        failed_rules = []

        for feature, rule in rules_config.items():
            val = row[feature]
            threshold = rule['thresh']
            operator = rule['op']

            passed = val > threshold if operator == '>' else val <= threshold

            if passed:
                score += 1
            else:
                failed_rules.append(feature)

        return score, failed_rules

    seq_df['Consensus_Score'], seq_df['Rules_Failed'] = zip(
        *seq_df.apply(calculate_consensus, axis=1)
    )

    # ==============================
    # 4. REPORT STRUCTURE
    # ==============================

    report_df = seq_df[['ID', 'Sequence', 'Consensus_Score']].copy()

    categories = ["10/10","9/10","8/10","7/10","6/10","5/10"]
    for cat in categories:
        target = int(cat.split('/')[0])
        report_df[cat] = report_df['Consensus_Score'].apply(
            lambda s: "✅" if s == target else ""
        )

    report_df['Failed_Descriptors'] = seq_df['Rules_Failed'].apply(
        lambda x: ", ".join(x) if x else "None"
    )

    ordered_features = [
        'VHSE5','charge','isoelectric_point','ProtFP4','SVGER5',
        'E4','T3','SVGER4','PRIN3','Z3'
    ]

    col_mapping = {}

    for feat in ordered_features:
        rule = rules_config[feat]
        col_name = f"{feat} ({rule['op']} {rule['thresh']})"
        report_df[col_name] = seq_df[feat]
        col_mapping[col_name] = feat

    report_df = report_df[report_df['Consensus_Score'] >= 5]
    report_df = report_df.sort_values(by='Consensus_Score', ascending=False)

    # ==============================
    # 5. STYLING
    # ==============================

    def highlight_failed(val, feature_key):
        rule = rules_config[feature_key]
        thresh = rule['thresh']
        op = rule['op']
        passed = val > thresh if op == '>' else val <= thresh
        if not passed:
            return 'color: red; font-weight: bold'
        return ''

    styled_report = report_df.style

    for col_name, feature_key in col_mapping.items():
        styled_report = styled_report.map(
            highlight_failed,
            feature_key=feature_key,
            subset=[col_name]
        )

    styled_report = styled_report.format(
        {col: "{:.3f}" for col in col_mapping.keys()}
    ).set_properties(**{'text-align': 'center'})

    # ==============================
    # DISPLAY
    # ==============================

    st.subheader("Consensus Report")
    st.dataframe(styled_report, use_container_width=True)

    # ==============================
    # DOWNLOAD
    # ==============================

    output = io.BytesIO()
    styled_report.to_excel(output, index=False, engine="openpyxl")
    output.seek(0)

    st.download_button(
        label="Download Excel Report",
        data=output,
        file_name=f"{name}_Consensus_Report.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )


# =========================
# About Us a tendina
# =========================
with st.expander("About Us"):
    st.markdown("""
    #### Original paper reference
    Davide Sestaioni, Simone Ventisette, Giulia Ciacci, Pasquale Palladino, Andrea Barucci, Maria Minunni, Simona Scarano (2026).  
    *Antibody-Free SPR Detection of Human Myoglobin in Serum by a Sandwich Configuration of Epitope-Imprinted Nanofilms and Nanoparticles*.  
    ACS Sensors 
    DOI: [doi link]
    Cite as:[]

    ---
    
    #### Research group
    Institute Name  
    University / CNR affiliation  

    Contact: email@domain.it
    """)
