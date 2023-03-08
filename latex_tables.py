def make_table(information,textfile_name,table_name):
    file = open ('data/'+textfile_name, 'w')
    file.write('\pgfplotstableread{')
    file.write('\n')   
    for i in range(len(information[0])):
        file.write(str(information[0][i])+ " "+ str(information[1][i]))
        file.write('\n')
    file.write("}{\ ".strip() + table_name+'}')

information =([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], [9.90578174668388, 9.31594925347879, 9.34793350184777, 9.35374706953787, 9.35290892923145, 9.35224179742688, 9.35229399518924, 9.35248577643269, 9.35250255803191, 9.35240297026792, 9.35236647664485, 9.35244534285036, 9.35250929440024, 9.35242492535370, 9.35229335172616])
textfile_name = 'Linear.txt'
table_name = 'LinearData'
make_table(information,textfile_name,table_name)