# teste_patentes.py
import requests
import json

url = "https://api.patentsview.org/patents/query"
payload = {
    "q": {"_text_any": ["title", "abstract"], "_search": "diabetes"},
    "f": ["patent_title", "patent_number"],
    "size": 5
}

try:
    r = requests.post(url, json=payload, timeout=10)
    dados = r.json()
    print("Status:", r.status_code)
    print("Total encontrado:", len(dados.get("patents", [])))
    for p in dados.get("patents", [])[:2]:
        print("Patente:", p.get("patent_title"))
except Exception as e:
    print("Erro:", e)