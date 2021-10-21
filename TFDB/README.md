# Transcription Factor (TF) Database

Database: ```sqlite3```

Tools: ```python flask```

Data type: ```json utf-8```

## information included in a piece of record
```
Accession：TEXT
Names：TEXT
Species：TEXT
Database：TEXT
Sequence：TEXT
Publications：TEXT
```

## HTTP Requests

### send requests： 

add data
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

query data
```
{
    "query" : "query",
    "data" : {
        "col_name" : key,
        "content" : value
    }
}
```

delete data
```
{
    "query" : "delete",
    "data" : {
        "col_name" : key,
        "content" : value
    }
}
```
### receive results：
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

