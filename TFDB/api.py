#coding=utf-8

from flask import Flask, request, jsonify
import os
import json
from database import RecordDAO
from flask_cors import CORS

trans_factor = RecordDAO("./trans_factor.db")

app = Flask(__name__)
CORS(app, supports_credentials=True, origin="http://sysu-software.com")

@app.route('/db', methods=['POST'])
def db():
    data = request.get_data()
    query = json.loads(data.decode())
    msg = {}
    if query['query'] == 'add':
        trans_factor.insert(query['data']['content'])
        msg['status'] = 0
        msg['message'] = 'succeed'
    elif query['query'] == 'query':
        if query['data']['col_name'] == "":
            raws = trans_factor.get_all()
        else:
            raws = trans_factor.get_condition(query['data']['col_name'], query['data']['content'])
        msg['data'] = {}
        msg['data']['col_name'] = trans_factor.col_name
        msg['data']['content'] = []
        for raw in raws:
            msg['data']['content'].append(raw[1:])
        msg['status'] = 0
        msg['message'] = 'succeed'
    elif query['query'] == 'delete':
        trans_factor.delete(query['data']['col_name'], query['data']['content'])
        msg['status'] = 0
        msg['message'] = 'succeed'
    else:
        msg['status'] = -1;
        msg['message'] = 'operation not allow'
    return jsonify(msg)
