import streamlit as st
import matplotlib.pyplot as plt
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image,
    PageBreak
)

from reportlab.lib.styles import getSampleStyleSheet
from Bio.Seq import Seq
from Bio.Restriction import *

from datetime import datetime

st.set_page_config(
    page_title="Restriction Enzyme Analysis Tool",
    page_icon="🧬",
    layout="wide"
)
# =====================================
# TITLE
# =====================================

st.title("🧬 Restriction Enzyme Analysis Tool")

# =====================================
# FASTA FILE UPLOAD
# =====================================

uploaded_file = st.file_uploader(
    "Upload FASTA File",
    type=["fasta", "fa", "txt"]
)

if uploaded_file is not None:

    content = uploaded_file.read().decode("utf-8")

    lines = content.splitlines()

    sequence = ""

    for line in lines:

        if not line.startswith(">"):
            sequence += line.strip()

    sequence = sequence.upper()

else:

    sequence = st.text_area(
        "Enter DNA Sequence"
    ).upper()



# =====================================
# ANALYZE BUTTON
# =====================================
enzyme_mode = st.selectbox(
    "Select Enzyme Database",
    [
        "Common Enzymes",
        "All Enzymes"
    ]
)
if st.button("Analyze"):

    st.session_state.run_analysis = True

if st.session_state.get("run_analysis", False):

    report_text = ""

    seq = Seq(sequence)

    if enzyme_mode == "Common Enzymes":
        analysis = Analysis(CommOnly, seq)
    else:
        analysis = Analysis(AllEnzymes, seq)

    results = analysis.full()

    available_enzymes = sorted(
        [
            str(enzyme)
            for enzyme, positions in results.items()
            if positions
        ]
    )

    if "selected_enzyme" not in st.session_state:
        st.session_state.selected_enzyme = available_enzymes[0]

    selected_enzyme = st.selectbox(
        "Select Enzyme",
        available_enzymes,
        index=available_enzymes.index(
            st.session_state.selected_enzyme
        )
    )

    st.session_state.selected_enzyme = selected_enzyme

    st.info(f"Using: {enzyme_mode}")

    left_col, right_col = st.columns(2)

    st.session_state.selected_enzyme = selected_enzyme
    
    total_sites = 0

    for enzyme, positions in results.items():

        total_sites += len(positions)

   
    st.info(
    f"Using: {enzyme_mode}"
)

    report_text += f"Sequence Length: {len(sequence)} bp\n\n"

    # =================================
    # ANALYSIS REPORT
    # =================================
    with left_col:
        st.subheader("Analysis Report")

        st.write(
            "Sequence Length:",
            len(sequence),
            "bp"
        )

        

        # =================================
        # FRAGMENT ANALYSIS
        # =================================

        st.subheader("Fragment Analysis")

        

        all_cuts = []

        for enzyme, positions in results.items():

            if str(enzyme) == selected_enzyme:

                all_cuts = positions.copy()
                break

        all_cuts.sort()

        fragments = []

        previous = 0

        for cut in all_cuts:

            fragments.append(cut - previous)

            previous = cut

        fragments.append(
            len(sequence) - previous
        )

        fragment_table = []
        MAX_FRAGMENTS = 200
        used_positions = []
      
        for i, fragment in enumerate(fragments[:MAX_FRAGMENTS]):

            cut_pos = ""

            if i < len(all_cuts):
                cut_pos = all_cuts[i]

            fragment_table.append(
                {
                    "Band": f"F{i+1}",
                    "Cut Position": cut_pos,
                    "Size (bp)": fragment
                }
            )



        st.table(fragment_table)
        report_text += (
            f"\nCut Positions: {all_cuts}\n"
        )

      

    with right_col:
        # =================================
        # RESTRICTION MAP
        # =================================

        st.subheader("Restriction Map")

        fig, ax = plt.subplots(figsize=(20,5))

        ax.plot(
            [0, len(sequence)],
            [0, 0],
            linewidth=3
        )

        for enzyme, positions in results.items():

            if str(enzyme) == selected_enzyme:

                for position in positions:

                    ax.vlines(
                        position,
                        0,
                        1
                    )

                    ax.text(
                        position,
                        1.05,
                        str(enzyme),
                        rotation=90,
                        fontsize=8
                    )

        ax.set_title("Restriction Map")

        ax.set_xlabel("DNA Position (bp)")

        ax.set_yticks([])

        fig.savefig(
            "restriction_map.png",
            bbox_inches="tight",
            dpi=300
        )

        st.pyplot(fig)

    

        # =================================
        # GEL ELECTROPHORESIS SIMULATION
        # =================================

        st.subheader("Gel Electrophoresis Simulation")



        gel_fig, gel_ax = plt.subplots(figsize=(7,10))

        gel_ax.set_facecolor("black")

        gel_ax.vlines(
            0,
            0,
            100,
            linewidth=25,
            color="white",
            alpha=0.15
        )

        import math

        band_positions = []

        for i, fragment in enumerate(fragments):

            position = 100 - (math.log10(fragment) * 25)

            if position < 5:
                position = 5

            # draw band
            gel_ax.hlines(
                position,
                -0.35,
                0.25,
                linewidth=5,
                color="white"
            )

            # group nearby bands
            found = False

            for band in band_positions:

                if abs(position - band["y"]) < 2:

                    band["labels"].append(f"F{i+1}")
                    found = True
                    break

            if not found:

                band_positions.append(
                    {
                        "y": position,
                        "labels": [f"F{i+1}"]
                    }
                )

        # draw grouped labels

        for band in band_positions:

            label_text = ", ".join(band["labels"])

            gel_ax.text(
                0.55,
                band["y"],
                label_text,
                color="yellow",
                fontsize=8,
                va="center",
                ha="left"
            )

        gel_ax.set_title("Simulated Agarose Gel")

        gel_ax.set_xticks([])
        gel_ax.set_yticks([])

        gel_ax.set_xlim(-1.0, 1.4)
        gel_ax.axvline(
            x=0.35,
            color="gray",
            alpha=0.3,
            linewidth=1
        )
        gel_fig.savefig(
            "gel_simulation.png",
            bbox_inches="tight",
            dpi=300
        )

        st.pyplot(gel_fig)

        
            # =================================
            # DOWNLOAD REPORT
            # =================================

        from reportlab.platypus import Table

        pdf_file = "restriction_report.pdf"

        doc = SimpleDocTemplate(pdf_file)

        styles = getSampleStyleSheet()

        elements = []

        elements.append(
            Paragraph(
                "Restriction Enzyme Analysis Report",
                styles["Title"]
            )
        )

        summary_table = Table([
            ["Parameter", "Value"],
            ["Analysis Date", datetime.now().strftime('%Y-%m-%d %H:%M')],
            ["Sequence Length", f"{len(sequence)} bp"],
            ["Selected Enzyme", selected_enzyme],
            ["Restriction Sites", str(len(all_cuts))],
            ["Fragments", str(len(fragments))]
        ])

        elements.append(summary_table)

        elements.append(Spacer(1,12))

        from datetime import datetime

        

        elements.append(Spacer(1,12))

        elements.append(
            Paragraph(
                f"Sequence Length: {len(sequence)} bp",
                styles["Normal"]
            )
        )

        elements.append(
            Paragraph(
                f"Selected Enzyme: {selected_enzyme}",
                styles["Normal"]
            )
        )

        elements.append(Spacer(1,12))

        elements.append(
            Paragraph(
                "Fragment Analysis",
                styles["Heading2"]
            )
        )

        table_data = [
            ["Fragment", "Size (bp)"]
        ]

        for i, fragment in enumerate(fragments):

            table_data.append(
                [
                    f"F{i+1}",
                    str(fragment)
                ]
            )

        fragment_table_pdf = Table(table_data)

        fragment_table_pdf.repeatRows = 1

        elements.append(fragment_table_pdf)
        
        elements.append(PageBreak())

        elements.append(
            Paragraph(
                "Restriction Map",
                styles["Heading2"]
            )
        )

        elements.append(
            Image(
                "restriction_map.png",
                width=500,
                height=150
            )
        )

        elements.append(Spacer(1,12))

        elements.append(
            Paragraph(
                "Gel Electrophoresis Simulation",
                styles["Heading2"]
            )
        )

        elements.append(
            Image(
                "gel_simulation.png",
                width=400,
                height=500
            )
        )

        elements.append(Spacer(1,20))

        elements.append(
            Paragraph(
                "Generated using Restriction Enzyme Analysis Tool v1.0",
                styles["Italic"]
            )
        )

        doc.build(elements)

        with open(
            pdf_file,
            "rb"
        ) as pdf:

            st.download_button(
                "📄 Download PDF Report",
                pdf,
                file_name="restriction_report.pdf",
                mime="application/pdf"
            )

            # =================================
            # DEBUG INFO
            # =================================

       

        
