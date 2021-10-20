#coding=utf-8

import sqlite3
import os

class RecordDAO:
    def __init__(self, db_path):
        self.db_path = db_path
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        self._create_table()
        self.col_name = self.get_col_name()

    def connect(self):
        return sqlite3.connect(self.db_path)

    def get_col_name(self):
        with self.connect() as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT * FROM protein"
            )
            col_name_list = [tuple[0] for tuple in cur.description if tuple[0] != "ID"]
        return col_name_list

    def _create_table(self):
        with self.connect() as conn:
            conn.execute(
                "CREATE TABLE IF NOT EXISTS protein ( \
                    ID INTEGER PRIMARY KEY AUTOINCREMENT, \
                    linker TEXT NOT NULL, \
                    长度 TEXT NOT NULL, \
                    所属全长域 INT NOT NULL, \
                    左侧保守域short_name TEXT NOT NULL, \
                    左侧保守域pssm_id TEXT NOT NULL, \
                    左侧保守域Accession TEXT NOT NULL, \
                    右侧保守域short_name TEXT NOT NULL, \
                    右侧保守域pssm_id TEXT NOT NULL, \
                    右侧保守域Accession TEXT NOT NULL,  \
                    溶液可及性SASA TEXT NOT NULL,  \
                    PBD文件编号PBD_ID TEXT NOT NULL  \
                    )"
            )

    def get_all(self):
        with self.connect() as conn:
            r = conn.execute(
                f"SELECT * FROM protein",
            ).fetchall()
        return r

    def condition_query(self, size, N, C):
        with self.connect() as conn:
            r = conn.execute(
                f"SELECT * FROM protein",
            ).fetchall()
        size = str(size)
        N = str(N)
        C = str(C)
        lsize = 0
        rsize = 0
        if(size == 'small'):
            lsize = 20
            rsize = 40
        elif(size == 'middle'):
            lsize = 40
            rsize = 60
        elif(size == 'large'):
            lsize = 40
            rsize = 60
        elif(size == 'very_large'):
            lsize = 80
            rsize = 100
        r_ans = []
        for _, node in enumerate(r):
            if(int(node[2]) >= lsize and int(node[2]) <= rsize):
                if(N == 'all' or str(node[1])[0] == N):
                    if(C == 'all' or str(node[1])[-1] == C):
                        r_ans.append(node)

        return r_ans

    def get_condition(self, key, value):
        key = str(key)
        value = "'" + str(value) + "'" if isinstance(value, str) else str(value)
        with self.connect() as conn:
            r = conn.execute(
                f"SELECT * FROM protein WHERE {key} = {value}",
            ).fetchall()
        return r

    def insert(self, datas):
        col = self.col_name
        col = ",".join(col)
        with self.connect() as conn:
            for data in datas:
                val = [f"'{i}'" if isinstance(i, str) else f"{i}" for i in data]
                val = ",".join(val)
                conn.execute(
                    f"INSERT INTO protein ({col}) VALUES ({val})"
                )

    def delete(self, key, value):
        key = str(key)
        value = "'" + str(value) + "'" if isinstance(value, str) else str(value)
        with self.connect() as conn:
            r = conn.execute(
                f"DELETE FROM protein where {key} = {value}",
            ).fetchall()
        return r