#R
#GJR 6/1/2018
#specifically queries 'RAD_linkage_matchPool' to find closest markers with linkage/chromosome information to a vector of scaffolds and positions

############## START ##################
queryDb.closeSnp<-function(scafs,pos,db='PomUrbanaGrant',Table='RAD_linkage_matchPool') {

    if ( length(scafs) != length(pos) ) {stop("lengths of scaffold and position vectors much match\n")}

    formQuery<-function(x) {
        scaf<-x[1]
        pos<-x[2]
        query<-paste("(
select   *
from     ",Table,"
where  scaffold = '",scaf,"'	 
and    position >= ",pos,"
order by position asc
limit 1
)
union
(
select   *
from     ",Table,"
where  scaffold = '",scaf,"'
and    position < ",pos,"
order by position desc
limit 1
)
order by abs(position - ",pos,")
limit 1;",sep="")
                                        #have to replace newlines which mysql doesn't like
        query<-gsub("\n", " ", query)
    }

    
    library(RMySQL)
    
    mydb = dbConnect(MySQL(), user='raglandlab', password='pomonella', dbname=db)
    j=1
    for (i in 1:length(scafs)) {
        query<-formQuery(c(scafs[i],pos[i]))
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
