#coding=utf-8

from flask import Flask, request, jsonify
from flask_cors import *
import os
import json
from database import RecordDAO

protein = RecordDAO("./linker_protein.db")

app = Flask(__name__)
CORS(app, supports_credentials=True)

@app.route('/db', methods=['POST'])
def db():
    data = request.get_data()
    query = json.loads(data.decode())
    msg = {}
    if query['query'] == 'add':
        protein.insert(query['data']['content'])
        msg['status'] = 0
        msg['message'] = 'succeed'
    elif query['query'] == 'Cquery':
        raws = protein.condition_query(query['data']['size'], query['data']['N'], query['data']['C'])
        msg['data'] = {}
        msg['data']['col_name'] = protein.col_name
        msg['data']['content'] = []
        for raw in raws:
            msg['data']['content'].append(raw[1:])
        msg['status'] = 0
        msg['message'] = 'succeed'
    elif query['query'] == 'query':
        if query['data']['col_name'] == "":
            raws = protein.get_all()
        else:
            raws = protein.get_condition(query['data']['col_name'], query['data']['content'])
        msg['data'] = {}
        msg['data']['col_name'] = protein.col_name
        msg['data']['content'] = []
        for raw in raws:
            msg['data']['content'].append(raw[1:])
        msg['status'] = 0
        msg['message'] = 'succeed'
    elif query['query'] == 'delete':
        protein.delete(query['data']['col_name'], query['data']['content'])
        msg['status'] = 0
        msg['message'] = 'succeed'
    else:
        msg['status'] = -1;
        msg['message'] = 'operation not allow'
    ans = jsonify(msg)
    ans.headers['Access-Control-Allow-Origin'] = '*'
    return ans