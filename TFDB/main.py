#coding=utf-8

from api import app

if __name__ == '__main__':
   app.config['JSON_AS_ASCII'] = False
   app.run(host='127.0.0.1', port='5122', debug=True)