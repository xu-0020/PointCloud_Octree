import psycopg2
import time
conn = psycopg2.connect(
    host="localhost",
    port="5433",
    database="hjingaa",
    user="postgres",
    password="postgres"
)
cur = conn.cursor()

file_path = "/export/project/hjingaa/PointCloud_Octree/result_hkm/Result_OR.txt"
sentence_list = []
count = 0
with open(file_path, "r") as file:
    lines = file.readlines()
    for line in lines:
        query_sentence = line.strip()
        tmp_minX = str(query_sentence.capitalize().split("(")[1].split(",")[0])
        tmp_maxX = str(query_sentence.capitalize().split(")")[0].split(" ")[-1])
        tmp_minY = str(query_sentence.capitalize().split("(")[2].split(",")[0])
        tmp_maxY = str(query_sentence.capitalize().split(")")[1].split(" ")[-1])
        tmp_minZ = str(query_sentence.capitalize().split("(")[3].split(",")[0])
        tmp_maxZ = str(query_sentence.capitalize().split(")")[2].split(" ")[-1])
        command = "select count(*) from hk_small where ST_3DIntersects(pos,Box3D(ST_SetSRID(ST_GeomFromEWKT('LINESTRING({tmp_minX} {tmp_minY} {tmp_minZ},{tmp_maxX} {tmp_maxY} {tmp_maxZ})'), 4326)));".format(tmp_minX=tmp_minX, tmp_maxX=tmp_maxX, tmp_minY=tmp_minY, tmp_maxY=tmp_maxY, tmp_minZ=tmp_minZ, tmp_maxZ=tmp_maxZ)
        # print(command)
        start_time = time.time()
        cur.execute(command)
        end_time = time.time()
        execution_time = end_time - start_time
        rows = cur.fetchall()
        for row in rows:
            print(row[0])
        sentence = "GISTGISTThe bound is (" + str(tmp_minX) + ", " + str(tmp_maxX) + "), (" + str(tmp_minY) + ", " + str(tmp_maxY) + "), (" + str(tmp_minZ) + ", " + str(tmp_maxZ) + ")." + " Time cost = " + str(execution_time) +" seconds." + "find " + str(row[0]) + " points.";
        sentence_list.append(sentence)
        print(sentence)
        # count += 1
        # if count == 3:
        #     break
cur.close()
conn.close()

with open("Result_tmp.txt", "w") as file:
    for item in sentence_list:
        file.write(item + "\n")