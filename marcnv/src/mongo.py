import pymongo



def get_mongo_database(uri: str, db_name: str):
    client = pymongo.MongoClient(uri)
    return client[db_name]
