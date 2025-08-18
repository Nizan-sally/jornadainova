# database.py

import os
import json
from sqlalchemy import create_engine, text
import hashlib
from dotenv import load_dotenv

load_dotenv()

# Monta a URL do banco 
DATABASE_URL = (
    f"postgresql://{os.getenv('SUPABASE_USER')}"
    f":{os.getenv('SUPABASE_PASS')}@"
    f"{os.getenv('SUPABASE_HOST')}:{os.getenv('SUPABASE_PORT')}/"
    f"{os.getenv('SUPABASE_DB')}"
)

def salvar_interacao(query, total_artigos=0, total_patentes=0):
    """
    Salva a busca com dados anonimizados
    """
    try:
        # echo=True debug
        engine = create_engine(DATABASE_URL, echo=False)
        fake_id = hashlib.sha256((os.getenv("USER_AGENT", "app") + query).encode()).hexdigest()[:10]
        sources_json = json.dumps({"artigos": total_artigos, "patentes": total_patentes})

        with engine.connect() as conn:
            conn.execute(
                text("""
                INSERT INTO user_interactions 
                (user_query, anonymized_id, results_count, data_sources)
                VALUES (:query, :user_id, :results, :sources)
                """),
                {
                    "query": query,
                    "user_id": fake_id,
                    "results": total_artigos + total_patentes,
                    "sources": sources_json
                }
            )
            conn.commit()
        print("Interação salva no banco")
    except Exception as e:
        print("Erro ao salvar no banco:", e)