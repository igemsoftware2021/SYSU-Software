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
                "SELECT * FROM trans_factor"
            )
            col_name_list = [tuple[0] for tuple in cur.description if tuple[0] != "ID"]
        return col_name_list

    def _create_table(self):
        with self.connect() as conn:
            conn.execute(
                "CREATE TABLE IF NOT EXISTS trans_factor ( \
                    ID INTEGER PRIMARY KEY AUTOINCREMENT, \
                    Accession TEXT NOT NULL, \
                    Names TEXT NOT NULL, \
                    Species TEXT NOT NULL, \
                    Database TEXT NOT NULL, \
                    Sequence TEXT NOT NULL, \
                    Publications TEXT NOT NULL \
                    )"
            )

    def get_all(self):
        with self.connect() as conn:
            r = conn.execute(
                f"SELECT * FROM trans_factor",
            ).fetchall()
        return r

    def get_condition(self, key, value):
        key = str(key)
        if key == 'Names':
            value = '%' + value + '%'
        value = "'" + str(value) + "'" if isinstance(value, str) else str(value)
        with self.connect() as conn:
            if key == 'Names':
                r = conn.execute(
                    f"SELECT * FROM trans_factor WHERE {key} LIKE {value}",
                ).fetchall()
            else:
                r = conn.execute(
                    f"SELECT * FROM trans_factor WHERE {key} = {value}",
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
                    f"INSERT INTO trans_factor ({col}) VALUES ({val})"
                )

    def delete(self, key, value):
        key = str(key)
        value = "'" + str(value) + "'" if isinstance(value, str) else str(value)
        with self.connect() as conn:
            r = conn.execute(
                f"DELETE FROM trans_factor where {key} = {value}",
            ).fetchall()
        return r