#R
#GJR 5/23/2018
#generic function to query tables in pomonella snp databases
#query string should be interpretable by sprintf to add 'id' as a character variable

############## START ##################
queryDb<-function(ids,queryString,db='PomUrbanaGrant') {


library(RMySQL)

mydb = dbConnect(MySQL(), user='raglandlab', password='pomonella', dbname=db)

j=1;
for (id in ids) {
    query<-sprintf(queryString,id)
    queryOb<-dbGetQuery(mydb, query)
    if (j==1) {geneList<-queryOb} else {
       geneList<-rbind(geneList,queryOb)
    }
    j=j+1
}
#terminate sql connections    
lapply(dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)
return(geneList)


}
############## END ##################








#example queries incorporating the id
#query<-sprintf("select * from annotation where snpId = %s", id)

#query<-sprintf("select * from feature_alias where loc ='%s'", id)


