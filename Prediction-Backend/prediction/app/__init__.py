from flask import Flask
import os

def after_request(resp):
    resp.headers['Access-Control-Allow-Origin'] = '*'
    return resp

ACTIVE_RESIDUAL_UPLOAD_FOLDER = '../userFiles/active/residual'
ACTIVE_RESIDUAL_ALLOWED_EXTENSIONS = {'pdb'}

ACTIVE_CAD_UPLOAD_FOLDER = '../userFiles/active/cad'
ACTIVE_CAD_ALLOWED_EXTENSIONS = {'pdb'}

STRUCTURE_UPLOAD_FOLDER = '../userFiles/structure'

SEQUENCE_UPLOAD_FOLDER = '../userFiles/sequence'
SEQUENCE_ALLOWED_EXTENSIONS = {'pdb', 'fa', 'fasta'}

TARGET_UPLOAD_FOLDER = '../userFiles/target'


app = Flask(__name__)
app.after_request(after_request)
app.config['ACTIVE_RESIDUAL_UPLOAD_FOLDER'] = ACTIVE_RESIDUAL_UPLOAD_FOLDER
app.config['ACTIVE_CAD_UPLOAD_FOLDER'] = ACTIVE_CAD_UPLOAD_FOLDER
app.config['STRUCTURE_UPLOAD_FOLDER'] = STRUCTURE_UPLOAD_FOLDER
app.config['SEQUENCE_UPLOAD_FOLDER'] = SEQUENCE_UPLOAD_FOLDER
app.config['TARGET_UPLOAD_FOLDER'] = TARGET_UPLOAD_FOLDER
app.config['SECRET_KEY'] = os.urandom(24)

from app import tasks
from app import views





