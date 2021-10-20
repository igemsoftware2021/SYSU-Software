# TF Database
这里是转录因子(TF)后端数据库的源码。

数据库：```sqlite3```

后端实现：```python flask```

数据类型：```json utf-8```

## 部署
```
pip install -r requirements.txt
python main.py
```

它运行在```localhost:5122/db```

## 数据表
```
Accession：TEXT
Names：TEXT
Species：TEXT
Database：TEXT
Sequence：TEXT
Publications：TEXT
```

## 数据格式
### 发送：
```
{
    "status" : 0 / -1,   
    "message" : "succeed",
    "data" : {
        "col_name" : [col1_name, col2_name, ...],
        "content" : [ [data1], [data2], ... ],
        ...
    }
}
```

### 接收： 

add
```
{
    "query" : "add",
    "data" : {
        "col_name" : [col1_name, col2_name, ...],
        "content" : [ [data1], [data2], ... ],
        ...
    }
}
```

query
```
{
    "query" : "query",
    "data" : {
        "col_name" : key,
        "content" : value
    }
}
```

delete,
```
{
    "query" : "delete",
    "data" : {
        "col_name" : key,
        "content" : value
    }
}
```

