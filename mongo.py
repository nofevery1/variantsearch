#!/usr/bin/python

from pymongo import MongoClient
import datetime
import re



#open mongo client
client = MongoClient("<add client address here>")
db = client['vsearch']
coll = db['variants']
updateFlag = None
docId = None


#is it already in db or not, consider date
def inDb(one,gene,rsid):
    searchDate = str(datetime.date.today())
    print one,gene,rsid
    if rsid is None:
        cursor = coll.find({"gene": gene, "protein": one})
    if one is None:
        cursor = coll.find({"gene": gene, "rsid": rsid})
    if rsid is not None and one is not None:
        cursor = coll.find({"gene": gene, "rsid": rsid, "protein": one})

    if cursor.count() == 0:
        print "no cursor"
        return False
    if cursor.count() > 0:
        if checkMissingFields(cursor,rsid,one) is True:
            global updateFlag
            updateFlag = "Update" #set this
            return False
        #if most recent document is more than a month old, then research pubmed, else return the articles from the document results
        for document in cursor.sort({"_id":-1}).limit(1):
            docDate = document['date']
        sDate = datetime.datetime.strptime(searchDate, "%Y-%m-%d").date()
        dDate = datetime.datetime.strptime(docDate, "%Y-%m-%d").date()
        delta = (sDate - dDate).days
        print delta
        #if the documents are older than one month, return to resarch pubmed and delete the deprecated documents before inserting new one
        if delta >= 30:
            deleted = coll.delete_many(
                { "$and": [
                    {"gene": gene},
                    {"$or": [
                        {"rsid": rsid},
                        {"protein": one}
                    ]}
                ]}
            )
            return False
        elif delta <30:
            return True



def checkMissingFields(cursor,rsid,one):
    #test if there are missing fields that can be entered, if there are then run the pubmed search again
    for document in cursor.sort("_id",-1).limit(1):
        print document
        global docId
        docId = document["_id"]
        rsidTest = str(document["rsid"])
        proteinTest = str(document["protein"])


    if re.match("None",rsidTest) and rsid is not None:
        updated = coll.update_one(
            {"_id": docId},
            {"$set": [
                {"rsid": rsid},
                ]}
        )
        return True
    elif re.match("None",proteinTest) and one is not None:
        updated = coll.update_one(
            {"_id": docId},
            {"$set": [
                {"protein": one},
                ]}
        )
        return True
    elif not re.match("None",rsidTest) and not re.match("None",proteinTest):
        return False



#insert document into database
def variantInput(one,gene,rsid,searchDict):
    searchDate = str(datetime.date.today())

    if rsid is None:
        rsid = "None"
    if one is None:
        one = "None"

    if updateFlag is None:
        result = coll.insert_one(
            {
                "gene": gene,
                "rsid": rsid,
                "protein": one,
                "date": searchDate,
                "paperCount": len(searchDict),
                "papers": searchDict
                }
            )
        return True
    elif re.match("Update",updateFlag):
        updated = coll.update_one(
            {"_id": docId},
            {"$set": [
                {"protein": one},
                {"date": searchDate},
                {"paperCount": len(searchDict)},
                {"papers": searchDict}
            ]}
        )
        return True


#if its already in the database, return the papers from the documents
def paperPull(one,gene,rsid):
    cursor = coll.find(
        { "$and": [
            {"gene": gene},
            {"$or": [
                {"rsid": rsid},
                {"protein": one}
            ]}
        ]}
    )

    docOut = dict()
    for document in cursor:
        docOut.update(document["papers"])
    return docOut
