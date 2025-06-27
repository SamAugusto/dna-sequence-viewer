import streamlit as st
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip
import io
import nest_asyncio
import random

# Required to avoid asyncio loop conflict in Streamlit
nest_asyncio.apply()

# Set NCBI email for legal API use
Entrez.email = "samaugusto121@gmail.com"  # ← Replace with your email

st.title("FASTQ/FASTA & NCBI Sequence Viewer")

# -------------------------
# Upload FASTQ Section
# -------------------------
uploaded_file = st.file_uploader("Upload a FASTQ file", type=["fastq", "fq", "gz",'fasta'])

def analyze_record(record):
    with st.expander(f"Read {record.id}", expanded=True):
        st.code(str(record.seq))
        data = Counter(record.seq)
        st.markdown("### Nucleotide Frequency")
        colors = ['#818D81', '#3C2A33', '#D3C9AC', '#E1702B']
        st.bar_chart(data,x_label = 'Nucleotides',y_label = 'Frequency',color = random.choice(colors))

        gc = 100 * gc_fraction(record.seq)
        st.write(f"**GC Content**: {gc:.2f}%")

        if "phred_quality" in record.letter_annotations:
            qualities = record.letter_annotations["phred_quality"]
            avg_quality = sum(qualities) / len(qualities)
            st.write(f"The average quality of Read **{record.id}** is: **{avg_quality:.2f}**")

# Process uploaded file
if uploaded_file is not None:
    st.success(f"Uploaded file: {uploaded_file.name}")

    if uploaded_file.name.endswith(".gz"):
        with gzip.open(uploaded_file, "rt") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fastq")):
                analyze_record(record)
                if i >= 50: break
    else:
        text_stream = io.TextIOWrapper(uploaded_file, encoding='utf-8')
        try:
            for i, record in enumerate(SeqIO.parse(text_stream, "fastq")):
                analyze_record(record)
                if i >= 50: break
        except Exception:
            text_stream.seek(0)  # Reset stream
            for i, record in enumerate(SeqIO.parse(text_stream, "fasta")):
                analyze_record(record)
                if i >= 50: break

# -------------------------
# NCBI Gene URL Search
# -------------------------
st.markdown("## Or Search by NCBI Gene URL")
search_url = st.text_input("Paste NCBI Gene URL here:")

def extract_gene_id_from_url(url):
    # Attempt to parse Gene ID from typical NCBI Gene URL
    import urllib.parse
    parsed = urllib.parse.urlparse(url)
    query = urllib.parse.parse_qs(parsed.query)
    term = query.get("Term") or query.get("id") or []
    return term[0] if term else None

def fetch_ncbi_fasta_by_gene_id(gene_id):
    try:
        # Search for linked nucleotide record(s)
        search_handle = Entrez.elink(dbfrom="gene", db="nucleotide", id=gene_id, linkname="gene_nuccore_refseqrna")
        link_result = Entrez.read(search_handle)
        search_handle.close()

        linksets = link_result[0]['LinkSetDb']
        if not linksets:
            return None

        nucleotide_id = linksets[0]['Link'][0]['Id']

        # Fetch the sequence from NCBI Nucleotide DB
        fetch_handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(fetch_handle, "fasta")
        fetch_handle.close()
        return seq_record
    except Exception as e:
        st.error(f"Error fetching from NCBI: {e}")
        return None

# Trigger fetch if a URL is entered
if search_url:
    gene_id = extract_gene_id_from_url(search_url)
    if gene_id:
        st.info(f"Fetching gene ID: {gene_id}")
        seq_record = fetch_ncbi_fasta_by_gene_id(gene_id)
        if seq_record:
            st.success(f"Successfully fetched: {seq_record.id}")
            analyze_record(seq_record)
        else:
            st.warning("No linked nucleotide sequence found for this gene.")
    else:
        st.error("Could not extract a valid gene ID from the URL.")
