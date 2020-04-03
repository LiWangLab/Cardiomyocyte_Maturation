

##README##

library(org.Mm.eg.db)
#keytypes(org.Hs.eg.db)
keytypes(org.Mm.eg.db)

load("demo.rdata")
############################################################################
####Datasets needs to be determined before running cross-match
####File1: the ENTREZID of differently expressed secreted proteins. Here in the demo you can find the cell_156_ligand_2
####File2: KEGG pathways' definitions and their related genes. Here in the demo you can find mmu_kegg_gene
####File3：a list of enriched pathways from multiple clusters. Here in the demo you can find enrich_sel
############################################################################

############################################################################
####Cross-match
############################################################################

####Find the KEGG pethways which the differently expressed secreted proteins are involved in.
cell_156_sec_id = AnnotationDbi::select(org.Mm.eg.db,keys= cell_156_ligand_2$gene,keytype="SYMBOL", columns = "ENTREZID", multiVals="first")  #  414 10%, 200 15%, 148 20%
cell_156_sec_kegg <- merge(unique(cell_156_sec_id),mmu_kegg_gene,by.x="ENTREZID",by.y="entrezgene")  #  578 10%, 494 15%, 334 20%  


####Then the pathways can be overlapped with the enrichment result

cross_gene_pathway_sel <- function (gene_term_set,gene_enrich_list,title,type="KEGG"){

	#gene_term_set is composed of two columns: ENTREZID and KEGG/GO ID
	#gene_enrich_list is list created by do_enrichment function.
    #type, KEGG, GO
	
	result<-data.frame()
	
	for (i in names(gene_enrich_list)){

		if(type=="KEGG"){
		
			tmp_merge<-merge(gene_term_set,gene_enrich_list[[i]][[1]],by.x="kegg_pathway",by.y="ID")
			
			if(nrow(tmp_merge)>0){
			
				new_title<-paste(c(title,i,"KG"),collapse="_")
				
				tmp_merge$cluster=new_title

				result<-rbind(result,tmp_merge)
			                     }
		
                        }
		
	}
	
	return(result)
}

cell_156_ligand_kegg_sel <- unique(cross_gene_pathway_sel(cell_156_sec_kegg, enrich_sel,"cell_156_sel"))  #   207 10%, 183 15%, 117 20% 

####Combine the overlapped pathways, secreted proteins and their respective cell clusters.
cell_156_ligand_kg <- cell_156_ligand_2[cell_156_ligand_2$gene %in% cell_156_ligand_kegg_sel$SYMBOL,] # 98 10%, 61 15%, 36 20%
cell_156_ligand_sel_kg <- merge(cell_156_ligand_kg, cell_156_ligand_kegg_sel, all=TRUE,by.x="gene", by.y="SYMBOL")[,c(1,7:10,16:20)]  # 604 10%, 221 15%, 127 20%
colnames(cell_156_ligand_sel_kg) <- c("Ligand","cell","ID","ENTREZID","definition","p.adjust","qvalue","GeneID","Count","Cluster")

cell_156_ligand_pathway <- cell_156_ligand_sel_kg  
cell_156_ligand_pathway$ligand_group <-paste(cell_156_ligand_pathway$cell, cell_156_ligand_pathway$Ligand,sep="_")  


#### To generate a network graph
library(network)
ligand_name <-sort(unique(cell_156_ligand_pathway$ligand_group))  # 98
term_name <-unique(cell_156_ligand_pathway$definition)   #29

cell_156_network <- matrix(0,nrow=length(ligand_name),ncol=length(term_name))
rownames(cell_156_network)<-ligand_name
colnames(cell_156_network)<-term_name

for (i in ligand_name){

	for (j in term_name){
	
		tmp_line<-subset(cell_156_ligand_pathway,ligand_group==i & definition == j)
		
		if(nrow(tmp_line)>0){
		
			cell_156_network[i,j]=1
		
		}
	
	}

}

cell_156_network<-network(cell_156_network,directed =TRUE)
cell_156_network_vertex<- network.vertex.names(cell_156_network)

cell_156_cell_type <- substr(ligand_name,1,5)
cell_156_network %v% "type" = ifelse(regexec("_",cell_156_network_vertex)!=-1, "ligand", "term")
cell_156_network %v% "cell_type" = c(cell_156_cell_type,rep("term",length(term_name)))
cell_156_network %v% "label_size" = ifelse(regexec("_",cell_156_network_vertex)!=-1, 3, 3)

library(GGally)
color_palette<- c("CM_2_" = "#FF6A6A", "CM_3_" = "#CD5555", "EC_1_" = "#63B8FF",
				  "FB_1_" = "#76EEC6", "FB_3_" = "#90EE90", "FB_4_" = "#43CD80","FB_5_"="#BCEE68",
				  "MP_3_" = "#FF8C00", "SMC_2" = "#FFFF00", "SMC_3" = "#FFD700", "term"="pink")

pdf("Interactionnetwork.pdf",height=10,width=12, family="ArialMT",useDingbats=FALSE)
p<-ggnet2(cell_156_network,size = "type",size.palette = c("term" = 18, "ligand" = 3),
         shape="type",shape.palette= c("term" = 18, "ligand" = 20),color="cell_type",
		 alpha = 1,edge.alpha = 1,palette=color_palette,
		 label.size="label_size",label= c(ligand_name,1:29))		 
#p<-font_theme(p)
print(p)
dev.off()


####Output the lists of pathways and secreted proteins 
cell_156_network_legend <- data.frame(node=1:29, pathway=ggnet2(cell_156_network,size = "type",
                           size.palette = c("term" = 6, "ligand" = 1),color="cell_type",
						   palette=color_palette,label.size="label_size",label=TRUE)$data [99:127,1])
write.xlsx(cell_156_network_legend, file="BR_cell_156_network_legend.xlsx", sheetName="correlation",row.names=FALSE)



