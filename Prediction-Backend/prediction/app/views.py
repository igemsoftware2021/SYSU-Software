from app import SEQUENCE_ALLOWED_EXTENSIONS, ACTIVE_RESIDUAL_ALLOWED_EXTENSIONS, ACTIVE_CAD_ALLOWED_EXTENSIONS, app
from app.tasks import get_candidate_residual, getCad, long_task
from flask import request, send_from_directory
from werkzeug.utils import secure_filename

import os
import sys
import shutil
import optipyzer as op

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser import PDBParser
import time


def active_residual_allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ACTIVE_RESIDUAL_ALLOWED_EXTENSIONS

def active_cad_allowed_file(src, tar):
    return '.' in src and \
        src.rsplit('.', 1)[1].lower() in ACTIVE_CAD_ALLOWED_EXTENSIONS and \
            '.' in tar and \
        tar.rsplit('.', 1)[1].lower() in ACTIVE_CAD_ALLOWED_EXTENSIONS

def sequence_allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in SEQUENCE_ALLOWED_EXTENSIONS

def fileId(filename):
    return filename[filename.rfind('/')+1 : filename.rfind('.') if filename.rfind('.') != -1 else len(filename)]

def parse_pdb(file_path, filename):
    fileSeq = SeqIO.parse(file_path, "pdb-atom")
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(filename[:filename.rfind('.')], file_path)
    first_chain_id = structure[0].get_list()[0].get_id()
    startPos = structure[0][first_chain_id].get_list()[0].get_id()[1]
    return fileSeq, startPos
    
def append_code(result_dict):
    if(result_dict.get('error')):
        return result_dict, 400
    else:
        return result_dict

@app.route('/hello')
def hello():
    return 'hello'

@app.route('/results/<pdbId>/<name>')
def download_file(pdbId, name):
    root_dir = os.path.dirname(os.getcwd())
    result_dir = os.path.join(os.path.join(root_dir, app.config['STRUCTURE_UPLOAD_FOLDER'][3:]), pdbId)
    return send_from_directory(result_dir, name, as_attachment=True)

@app.route('/active/residual', methods=['POST'])
def active_residual():
    result_dict = {}
    # check if the post request has the file part
    if 'file' not in request.files:
        result_dict['error'] = 'No file part'
        return append_code(result_dict)
    file = request.files['file']
    if file.filename =='':
        result_dict['error'] = 'No selected file'
        return append_code(result_dict)
    if file and active_residual_allowed_file(file.filename):
        filename = secure_filename(file.filename)
        pdbId = filename[:filename.rfind('.')] + '_' + time.strftime('%Y-%m-%d-%H-%M-%S',time.localtime(time.time()))
        file.save(os.path.join(app.config['ACTIVE_RESIDUAL_UPLOAD_FOLDER'], filename))
        result = {}
        try:
            result['active_list'] = get_candidate_residual(os.path.join(app.config['ACTIVE_RESIDUAL_UPLOAD_FOLDER'], filename))
        except:
            error_title = str( sys.exc_info()[0] )
            error_detail = str( sys.exc_info() )
            result_dict['error'] = error_title + ": " + error_detail
            return append_code(result_dict)
        dest_dir = os.path.join(app.config['ACTIVE_CAD_UPLOAD_FOLDER'], pdbId)
        if not os.path.exists(dest_dir):
            os.mkdir(dest_dir)
        file_path = os.path.join(dest_dir, filename)
        shutil.move(os.path.join(app.config['ACTIVE_RESIDUAL_UPLOAD_FOLDER'], filename), file_path)
        result['srcPath'] = os.path.join(pdbId, filename)
        result_dict['result'] = result
        fileSeq, startPos = parse_pdb(file_path, filename)
        resultSeq = ''
        fileSeq = list(fileSeq)
        
        for record in fileSeq:
            resultSeq += str(record.seq)
            break
        
        result['seq'] = resultSeq
        result['startPos'] = startPos
    else:
        result_dict['error'] = 'The extension of file is not allowed.'
    return append_code(result_dict)
    
@app.route('/active/cad', methods=['POST'])
def active_cad():
    result_dict = {}
    srcPath = request.form.get('srcPath')
    tarPath = request.form.get('tarPath')
    print(srcPath)
    print(tarPath)
    # srcPath = srcPath.replace('%2F', '/')
    # tarPath = tarPath.replace('%2F', '/')

    if not srcPath or (not tarPath and 'tarPDB' not in request.files):
        result_dict['error'] = 'Not enough input.'
        return append_code(result_dict)
    
    try:
        srcPDBname = srcPath[srcPath.rfind('/')+1:]
        srcPath = srcPath[:srcPath.rfind('/')]
        pdbPath = os.path.join(app.config['ACTIVE_CAD_UPLOAD_FOLDER'], srcPath)
        tarPDBname = 'fusion_'
        
        if not tarPath:
            tarPDB = request.files['tarPDB']
            tarPDBname += secure_filename(tarPDB.filename)
            tarPDB.save(os.path.join(pdbPath, tarPDBname))
        else:
            tarPDBname += tarPath[tarPath.rfind('/')+1:]
            shutil.copyfile(os.path.join(app.config['STRUCTURE_UPLOAD_FOLDER'], tarPath), os.path.join(pdbPath, tarPDBname))
        
        active_list = request.form.getlist('active_list[]')
        thresh = request.form.get('thresh')
        if thresh:
            thresh = int(thresh)
 
        cad = getCad(pdbPath, srcPDBname, tarPDBname, active_list, thresh)
        if cad == -1:
            result_dict['error'] = "The source pdb file and target pdb file mismatch."
            return append_code(result_dict)
        elif cad == -2:
            result_dict['error'] = "Not enough valid active residual to calculate cad."
            return append_code(result_dict)
        if not tarPath:
            cad['tarPath'] = os.path.join(srcPath, tarPDBname)
        result_dict['result'] = cad
    except:
        error_title = str( sys.exc_info()[0] )
        error_detail = str( sys.exc_info() )
        result_dict['error'] = error_title + ": " + error_detail
        # os.system('rm %s -r'%pdbPath)
        return append_code(result_dict)
    # os.system('rm %s -r'%pdbPath)
    return append_code(result_dict)

@app.route('/active/<path>/<name>')
def get_active_file(path, name):
    root_dir = os.path.dirname(os.getcwd())
    active_dir = os.path.join(root_dir, app.config['ACTIVE_CAD_UPLOAD_FOLDER'][3:])
    return send_from_directory(os.path.join(active_dir, path), name, as_attachment=True)

@app.route('/active/quit', methods=['POST'])
def avtive_quit():
    srcPath = request.form.get('srcPath')
    srcPath = os.path.join(app.config['ACTIVE_CAD_UPLOAD_FOLDER'], srcPath[:srcPath.rfind('/')])
    if os.path.exists(srcPath):
        os.system('rm %s -r'%srcPath)
    return 'success'

@app.route('/structure', methods=['POST'])
def structure_predict():
    result_dict = {}
    linker = request.form.get('linker')
    lightSeq = request.form.get('lightSeq')
    tarSeq = request.form.get('tarSeq')
    id = request.form.get('id')

    if not linker or not lightSeq or not id or not tarSeq:
        result_dict['error'] = 'Not enough input.'
        return append_code(result_dict)

    pdbPath = os.path.join(app.config['STRUCTURE_UPLOAD_FOLDER'], id)
    if not os.path.exists(pdbPath):
        os.mkdir(pdbPath)

    resultSeq = ''
    resultSeq += lightSeq
    resultSeq += linker
    resultSeq += tarSeq
    
    max_len = 600
    resultSeq = ''.join(list(filter(lambda ch : not (ch == 'X'), resultSeq)))
    if len(resultSeq) > max_len or len(resultSeq) < 27:
        result_dict['error'] = 'Proper sequence length: 26 < length < 601.'
        return append_code(result_dict)

    resultRec = SeqRecord(
        Seq(
            resultSeq
        ),
        id=id,
        description='Light-activated protein Linker Target protein',
    )

    SeqIO.write([resultRec], os.path.join(pdbPath, id + ".fa"), "fasta")
    task_id = long_task.delay(pdbPath, id + ".fa", id).task_id
    result_dict['task_id'] = task_id
    return append_code(result_dict)

@app.route('/get/<task_id>')
def get(task_id):
  result = long_task.AsyncResult(task_id)
  status = result.status
  info = result.info
#   result.forget()
  result_dict = {}
  result_dict['status'] = status
  result_dict['info'] = info
  return append_code(result_dict)

@app.route('/sequence', methods=['POST'])
def sequence():
    result_dict = {}
    # check if the post request has the file part
    if 'file' not in request.files:
        result_dict['error'] = 'No file part'
        return append_code(result_dict)
    file = request.files['file']
    if file.filename =='':
        result_dict['error'] = 'No selected file'
        return append_code(result_dict)
    if file and sequence_allowed_file(file.filename):
        filename = secure_filename(file.filename)
        file_path = os.path.join(app.config['SEQUENCE_UPLOAD_FOLDER'], filename)
        file.save(file_path)
        try:
            ext = filename[filename.rfind('.') + 1:]
            if ext.startswith('fa'):
                fileSeq = SeqIO.parse(file_path, "fasta")
                startPos = 1
            else:
                fileSeq, startPos = parse_pdb(file_path, filename)

            resultSeq = ''
            fileSeq = list(fileSeq)
            # if len(fileSeq) > 1:
            #     result_dict['warning'] = 'Not support multichain. Only provide the first chain.'
            for record in fileSeq:
                resultSeq += str(record.seq)
                break
            result = {}
            result['seq'] = resultSeq
            result['id'] = filename[:filename.rfind('.')]
            result['startPos'] = startPos
            
        except:
            error_title = str( sys.exc_info()[0] )
            error_detail = str( sys.exc_info() )
            result_dict['error'] = error_title + ": " + error_detail
            return append_code(result_dict)
        save = request.form.get('save')
        if save and save == 'target':
            shutil.move(file_path, os.path.join(app.config['TARGET_UPLOAD_FOLDER'], filename))
            result['tarPath'] = filename
        else:
            os.remove(file_path)
    else:
        result_dict['error'] = 'The extension of file is not allowed.'
    result_dict['result'] = result
    return append_code(result_dict)

@app.route('/target/<name>')
def get_target_file(name):
    root_dir = os.path.dirname(os.getcwd())
    tar_dir = os.path.join(root_dir, app.config['TARGET_UPLOAD_FOLDER'][3:])
    return send_from_directory(tar_dir, name, as_attachment=True)

@app.route('/optimize', methods=['POST'])
def optimize():
    result_dict = {}
    optipyzer = op.api()
    species_list = request.form.getlist('species_list[]')
    weights_list = request.form.getlist('weights_list[]')
    seq = request.form.get('seq')
    weights_list = [int(x) for x in weights_list]
    
    if not seq or len(species_list) == 0 or len(weights_list) == 0:
        result_dict['error'] = 'Not enough input.'
        return append_code(result_dict)

    if seq.find('X') != -1:
        result_dict['error'] = 'The protein sequence contains a non-residue X.'
        return append_code(result_dict)

    try:
        print(seq)
        print(weights_list)
        print(species_list)
        org_list = []
        codon_usage_list = []
        for specie in species_list:
            results = optipyzer.search(name=specie)
            org = results[0]
            org_list.append(org)
            codon_usage_list.append(optipyzer.pull_codons(org))

        optimized = optipyzer.optimize(seq,org_list=org_list,weights=weights_list, seq_type='protein')
        result = {}
        result['peptide_seq'] = optimized.peptide_seq
        result['optimmized_sd'] = optimized.optimmized_sd
        result['optimmized_ad'] = optimized.optimmized_ad
        
    except:
        error_title = str( sys.exc_info()[0] )
        error_detail = str( sys.exc_info() )
        result_dict['error'] = error_title + ": " + error_detail
        return append_code(result_dict)

    
    result_dict['result'] = result
    return append_code(result_dict)

@app.route('/test1')
def test():
    return 'test'