# teste_pubmed.py
from Bio import Entrez
import ssl

# Ignora certificado (só para teste)
ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "teste@exemplo.com"

try:
    handle = Entrez.esearch(db="pubmed", term="diabetes", retmax=5)
    result = Entrez.read(handle)
    handle.close()
    print(" Sucesso! Artigos encontrados:", result["IdList"])
except Exception as e:
    print("Falha no PubMed:", e)