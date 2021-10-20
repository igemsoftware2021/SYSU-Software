# -*- coding: utf-8 -*-
import os
from flask import Flask,  send_from_directory, request, jsonify, Response
from werkzeug.utils import secure_filename
import json
from EnzymeActiveSites import calculateActiveParam
from lightVersusRate import lightVersusRateOfReaction
from flask_cors import CORS, cross_origin


app = Flask(__name__)
CORS(app, supports_credentials=True, origin="http://sysu-software.com")

UPLOAD_PDB = 'upload_pdb'
UPLOAD_JSON = 'upload_json'
app.config['UPLOAD_PDB'] = UPLOAD_PDB  
app.config['UPLOAD_JSON'] = UPLOAD_JSON  
basedir = os.path.abspath(os.path.dirname(__file__))  
ALLOWED_EXTENSIONS = set(['pdb'])  

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def gen_filename(filename, pdb = True):
    i = 1
    type = app.config['UPLOAD_PDB'] if pdb else app.config['UPLOAD_JSON']
    sub = "pdb" if pdb else "json"
    filename = filename.split('.')[1] + '.' + sub
    newfilename = filename
    while os.path.exists(os.path.join(type, newfilename)):
        newfilename = f"{i}_{filename}"
        i += 1
    return newfilename

@app.route('/api/mathmodel/calcactiveparam', methods=['POST'], strict_slashes=False)
def api_calcactiveparam():
    file_dir = os.path.join(basedir, app.config['UPLOAD_PDB'])  
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)  
    pdbs = []
    pdbs.append(request.files['pdb1'])
    pdbs.append(request.files['pdb2'])
    pdbs.append(request.files['pdb3'])
    pdbs.append(request.files['pdb12'])
    pdbs.append(request.files['pdb23'])
    pdbname = []
    for pdb in pdbs:
        filename = secure_filename(pdb.filename)
        if allowed_file(filename):  
            filename = gen_filename(filename)
            filename = os.path.join(file_dir, filename)
            pdbname.append(filename)
            pdb.save(filename)  
        else:
            msg = {}
            msg['message'] = "非法文件后缀名！"
            return jsonify(msg), 400
    
    sname = []
    sname.append(request.form['sjson1'])
    sname.append(request.form['sjson2'])
    sname.append(request.form['sjson3'])

    try:
        sIntersection1, alpha1, r1, sIntersection2, alpha2, r2, d12 = calculateActiveParam(pdbname[0], pdbname[1], pdbname[3], sname[0], sname[1])
    except Exception as e:
        s = str(e)
        return jsonify({"err" : s}), 500

    datas = {}
    datas['result1'] = {}
    datas['result1']['sIntersection1'] = sIntersection1
    datas['result1']['sIntersection2'] = sIntersection2
    datas['result1']['d12'] = d12

    try:
        sIntersection2, alpha2, r2, sIntersection3, alpha3, r3, d23 = calculateActiveParam(pdbname[1], pdbname[2], pdbname[4], sname[1], sname[2])
    except Exception as e:
        s = str(e)
        return jsonify({"err" : s}), 500

    datas['result2'] = {}
    datas['result2']['sIntersection2'] = sIntersection2
    datas['result2']['sIntersection3'] = sIntersection3
    datas['result2']['d23'] = d23
    return jsonify(datas), 200
            


@app.route('/api/mathmodel/lightversusrateofreaction', methods=['POST'], strict_slashes=False)
def api_lightversusrateofreaction():
    datas = request.get_json()
    a1 = eval(datas['a1'])
    d1 = eval(datas['d1'])
    k1 = eval(datas['k1'])
    a2 = eval(datas['a2'])
    d2 = eval(datas['d2'])
    k2 = eval(datas['k2'])
    a3 = eval(datas['a3'])
    d3 = eval(datas['d3'])
    k3 = eval(datas['k3'])
    s = eval(datas['s'])
    etol = eval(datas['etol'])
    lkFuncStr = str(datas['lkFuncStr'])
    k1m = eval(datas['k1m'])
    k2m = eval(datas['k2m'])
    A1 = float(datas['A1'])
    A2 = float(datas['A2'])
    l12 = float(datas['l12'])
    l23 = float(datas['l23'])
    D1 = eval(datas['D1'])
    D2 = eval(datas['D2'])
    try:
        fittedFunc, graphFilename, fittedValue = lightVersusRateOfReaction(a1, d1, k1, a2, d2, a3, d3, k2, k3, s, etol, lkFuncStr, k1m, k2m, A1, A2, l12, l23, D1, D2)
    except Exception as e:
        s = str(e)
        return jsonify({"err" : s}), 500
    msg = {}
    msg['pic'] = graphFilename
    keys = [str(x) for x in range(len(fittedFunc))]
    keys = keys[::-1]
    list_json = dict(zip(keys, fittedFunc))
    msg['result'] = list_json
    msg['resultarr'] = fittedValue.tolist()
    return jsonify(msg), 200

@app.route("/api/mathmodel/img/<filename>")
def downloader(filename="README.md"):
    dirpath = app.root_path
    # dirpath = os.path.join(app.root_path, 'upload')
    return send_from_directory(dirpath, filename), 200


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8101, debug=True)
