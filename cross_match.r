
############################################################################
####Datasets needs to be determined before running cross-match
####File1: the ENTREZID of differently expressed secreted proteins. Here in the demo you can find the cell_156_ligand_2
####File2: KEGG pathways' definitions and their related genes. Here in the demo you can find mmu_kegg_gene
####File2 is created following the steps below:
    #1. The KEGG pathways associated with genes were retrieved by the program kg (https://github.com/endrebak/kg) 
        kg -s mmu -d >mmu_kegg_path2gene

       #-s SPEC --species=SPEC  name of species (examples: hsa, mmu, rno...) 
       #-d --definitions        add KEGG pathway definitions to the output 

    #2. add mmu before pathway ID
       awk -F "\t" '{if ($1 !~/pathway/){ $1="mmu"$1;} print $1"\t"$2"\t"$3}' mmu_kegg_path2gene >mmu_kegg_path2gene_full
       
	#3. read the file into R 
	   mmu_kegg_gene<-read.delim("mmu_kegg_path2gene_full",header=T,as.is=T)
	   
####File3：a list of enriched pathways from multiple clusters. Here in the demo you can find enrich_sel
############################################################################


############################################################################
####Cross-match
############################################################################

####Open a R 3.5 session
library(xlsx)
library(AnnotationDbi) 
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)

load("demo.rdata")

####Find the KEGG pathways which the differently expressed secreted proteins are involved in.
sec_id = AnnotationDbi::select(org.Mm.eg.db,keys= cell_156_ligand_2$gene,keytype="SYMBOL", columns = "ENTREZID", multiVals="first")  
sec_kegg <- merge(unique(sec_id),mmu_kegg_gene,by.x="ENTREZID",by.y="entrezgene")  


####Then the pathways can be overlapped with the enrichment result.

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

ligand_kegg_sel <- unique(cross_gene_pathway_sel(sec_kegg, enrich_sel,"cell_156_sel"))  

####Combine the overlapped pathways, secreted proteins and their respective cell clusters.
ligand_kg <- cell_156_ligand_2[cell_156_ligand_2$gene %in% ligand_kegg_sel$SYMBOL,] 
ligand_sel_kg <- merge(ligand_kg, ligand_kegg_sel, all=TRUE,by.x="gene", by.y="SYMBOL")[,c(1,7:10,16:20)]  
colnames(ligand_sel_kg) <- c("Ligand","cell","ID","ENTREZID","definition","p.adjust","qvalue","GeneID","Count","Cluster")

ligand_pathway <- ligand_sel_kg  
ligand_pathway$ligand_group <-paste(ligand_pathway$cell, ligand_pathway$Ligand,sep="_")  


####To generate a network graph. Please NOTE: the spatial structure of the network can be different from each time the network is produced. 
library(network)
ligand_name <-sort(unique(ligand_pathway$ligand_group))  
term_name <-unique(ligand_pathway$definition)   

cell_network <- matrix(0,nrow=length(ligand_name),ncol=length(term_name))
rownames(cell_network)<-ligand_name
colnames(cell_network)<-term_name

for (i in ligand_name){

	for (j in term_name){
	
		tmp_line<-subset(ligand_pathway,ligand_group==i & definition == j)
		
		if(nrow(tmp_line)>0){
		
			cell_network[i,j]=1
		
		}
	
	}

}

cell_network<-network(cell_network,directed =TRUE)
cell_network_vertex<- network.vertex.names(cell_network)

cell_type <- substr(ligand_name,1,5)
cell_network %v% "type" = ifelse(regexec("_",cell_network_vertex)!=-1, "ligand", "term")
cell_network %v% "cell_type" = c(cell_type,rep("term",length(term_name)))
cell_network %v% "label_size" = ifelse(regexec("_",cell_network_vertex)!=-1, 3, 3)

library(GGally)

##Set the color
color_palette<- c("CM_2_" = "#FF6A6A", "CM_3_" = "#CD5555", "EC_1_" = "#63B8FF",
				  "FB_1_" = "#76EEC6", "FB_3_" = "#90EE90", "FB_4_" = "#43CD80","FB_5_"="#BCEE68",
				  "MP_3_" = "#FF8C00", "SMC_2" = "#FFFF00", "SMC_3" = "#FFD700", "term"="pink")

##The label=c() should be set
pdf("Interactionnetwork.pdf",height=10,width=12, family="ArialMT",useDingbats=FALSE)
p<-ggnet2(cell_network,size = "type",size.palette = c("term" = 18, "ligand" = 3),
         shape="type",shape.palette= c("term" = 18, "ligand" = 20),color="cell_type",
		 alpha = 1,edge.alpha = 1,palette=color_palette,
		 label.size="label_size",label= c(ligand_name,1:29))		 
print(p)
dev.off()


####Output the lists of pathways and secreted proteins 
cell_network_legend <- data.frame(node=1:29, pathway=ggnet2(cell_network,size = "type",
                           size.palette = c("term" = 6, "ligand" = 1),color="cell_type",
						   palette=color_palette,label.size="label_size",label=TRUE)$data [99:127,1])
write.xlsx(cell_network_legend, file="network_legend.xlsx", sheetName="correlation",row.names=FALSE)


############################################################################
####Session Info
############################################################################

R version 3.5.0 (2018-04-23)

org.Mm.eg.db_3.6.0   AnnotationDbi_1.42.1 IRanges_2.14.12
S4Vectors_0.18.3     Biobase_2.40.0       BiocGenerics_0.26.0
xlsx_0.6.1           Rcpp_1.0.1           xlsxjars_0.6.1  
digest_0.6.25        DBI_1.0.0            RSQLite_2.1.1   
blob_1.1.1           bit64_0.9-7          bit_1.1-14
compiler_3.5.0       pkgconfig_2.0.3      rJava_0.9-10    
memoise_1.1.0

