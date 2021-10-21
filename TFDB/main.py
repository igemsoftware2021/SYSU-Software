#coding=utf-8

from api import app

if __name__ == '__main__':
   app.config['JSON_AS_ASCII'] = False
   app.run(host='0.0.0.0', port='5122', debug=True)