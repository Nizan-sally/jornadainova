# data_fetcher.py
import ssl, os, random
import pandas as pd
from Bio import Entrez
from dotenv import load_dotenv
from googletrans import Translator
import streamlit as st
from keybert import KeyBERT

ssl._create_default_https_context = ssl._create_unverified_context
load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL")

# Tradutor síncrono
translator = Translator()
kw_model = KeyBERT()

STOPWORDS = ["de", "do", "da", "em", "para", "com", "e", "o", "a", "os", "as"]
STOPWORDS_SAUDE = ["paciente","criança","homem","mulher","idoso","adulto",
                   "caso","grupo","estudo","dados","ano","resultado"]

# construir query
def construir_query(tema, foco=None, data_inicio=None, data_fim=None):
    palavras = [p for p in tema.split() if p.lower() not in STOPWORDS] or [tema]
    query = " OR ".join([f"{p}[Title/Abstract]" for p in palavras])
    focos = {
        "Tratamento": "(therapy OR treatment OR intervention)",
        "Diagnóstico": "(diagnosis OR detection OR screening)",
        "Prevenção": "(prevention OR prophylaxis OR vaccine)"
    }
    if foco in focos:
        query += f" AND {focos[foco]}[Title/Abstract]"
    if data_inicio or data_fim:
        query += f" AND ({data_inicio or 1800}[Date - Publication] : {data_fim or 3000}[Date - Publication])"
    return query

def extrair_info(article):
    ano = "N/A"
    if "ArticleDate" in article and len(article["ArticleDate"]) > 0:
        ano = article["ArticleDate"][0].get("Year", "N/A")
    elif "Journal" in article:
        ano = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {}).get("Year", "N/A")
    autores = ", ".join([f"{a.get('LastName','N/A')} {a.get('Initials','')}" for a in article.get("AuthorList", [])[:2]]) or "N/A"
    return ano, autores

@st.cache_data(show_spinner="Buscando artigos...", allow_output_mutation=True)
def buscar_artigos(tema, max_itens=15, foco=None, data_inicio=None, data_fim=None):
    try:
        query = construir_query(tema, foco, data_inicio, data_fim)
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_itens)
        ids = Entrez.read(handle).get("IdList", [])
        handle.close()
        if not ids:
            return pd.DataFrame(columns=["Título","Autores","Ano"])

        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        artigos = []
        for rec in records.get("PubmedArticle", []):
            art = rec.get("MedlineCitation", {}).get("Article", {})
            title = art.get("ArticleTitle", "Sem título")
            try:
                # Tradução síncrona
                title_pt = translator.translate(title, src="auto", dest="pt").text if title else title
            except:
                title_pt = title
            ano, autores = extrair_info(art)
            artigos.append({"Título": title_pt, "Autores": autores, "Ano": ano})

        return pd.DataFrame(artigos, columns=["Título","Autores","Ano"])
    except Exception as e:
        print("[Erro PubMed]:", e)
        return pd.DataFrame(columns=["Título","Autores","Ano"])
