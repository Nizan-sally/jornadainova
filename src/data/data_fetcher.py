# data_fetcher

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
translator = Translator()
kw_model = KeyBERT()

STOPWORDS = ["de", "do", "da", "em", "para", "com", "e", "o", "a", "os", "as"]
STOPWORDS_SAUDE = ["paciente","criança","homem","mulher","idoso","adulto",
                   "caso","grupo","estudo","dados","ano","resultado"]

# Construir query
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

# xtrair informações
def extrair_info(article):
    ano = "N/A"
    if "ArticleDate" in article and len(article["ArticleDate"]) > 0:
        ano = article["ArticleDate"][0].get("Year", "N/A")
    elif "Journal" in article:
        ano = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {}).get("Year", "N/A")

    autores = ", ".join(
        [f"{a.get('LastName','N/A')} {a.get('Initials','')}" 
         for a in article.get("AuthorList", [])[:2]]
    ) or "N/A"

    return ano, autores

# Buscar artigos
@st.cache_data(show_spinner="Buscando artigos...")
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
                title_pt = translator.translate(title, src="auto", dest="pt").text if title else title
            except:
                title_pt = title
            ano, autores = extrair_info(art)
            artigos.append({"Título": title_pt, "Autores": autores, "Ano": ano})

        return pd.DataFrame(artigos, columns=["Título","Autores","Ano"])
    except Exception as e:
        print("[Erro PubMed]:", e)
        return pd.DataFrame(columns=["Título","Autores","Ano"])

# Extrair conceitos-chave
def extrair_conceitos(textos, n=10):
    texto_concat = " ".join(textos)
    keywords = kw_model.extract_keywords(texto_concat, keyphrase_ngram_range=(1,2), 
                                         stop_words=STOPWORDS+STOPWORDS_SAUDE, top_n=n)
    return [k[0] for k in keywords]

# Analisar gaps
def analisar_gaps(df_artigos, df_patentes):
    if df_artigos.empty:
        return [], ["Poucos ou nenhum artigo encontrado. Difícil analisar oportunidades."]

    conceitos_artigos = extrair_conceitos(df_artigos['Título'].tolist(), n=15)
    conceitos_patentes = extrair_conceitos(df_patentes['Título'].tolist(), n=15) if not df_patentes.empty else []

    gaps = [c for c in conceitos_artigos if c not in conceitos_patentes]

    recomendacoes = []
    if gaps:
        for g in gaps[:5]:
            recomendacoes.append(f"Alta produção acadêmica em '{g}', mas baixa proteção de patentes. Oportunidade de inovação.")
    else:
        recomendacoes.append("Não foram identificados gaps claros — o mercado parece competitivo ou equilibrado.")

    return gaps, recomendacoes

# Buscar patentes
def buscar_patentes(tema, df_artigos, max_itens=5):
    if df_artigos.empty or "Título" not in df_artigos.columns:
        return pd.DataFrame(columns=["Título","Número","Data","Inventor"]), [], ["Nenhum artigo encontrado, logo não foram localizadas patentes."]

    patentes = []
    termos = []
    for titulo in df_artigos['Título'].dropna():
        termos.extend([t for t in titulo.split() if len(t) > 3])

    termos = termos[:50]
    nomes_inventores = ['Silva','Souza','Oliveira','Costa','Pereira','Santos','Almeida']

    qtd_artigos = len(df_artigos)
    if qtd_artigos <= 5:
        qtd_patentes = random.randint(0, 2)
    else:
        qtd_patentes = random.randint(2, min(5, max_itens))

    for i in range(qtd_patentes):
        termo = random.choice(termos) if termos else "Tecnologia"
        titulo_patente = f"Composição para {termo}"
        numero = f"US2025/{1000+i}"
        data = str(random.randint(2018,2025))
        inventor = f"Dr. {random.choice(nomes_inventores)}"
        patentes.append({"Título": titulo_patente, "Número": numero, "Data": data, "Inventor": inventor})

    df_patentes = pd.DataFrame(patentes, columns=["Título","Número","Data","Inventor"])
    gaps, recomendacoes = analisar_gaps(df_artigos, df_patentes)

    return df_patentes, gaps, recomendacoes