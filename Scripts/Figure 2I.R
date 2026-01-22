library(openxlsx)
library(Loafer) ## install_github("15652939484/Loafer")
library(limma)
library(pheatmap)
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)
library(magrittr)
library(igraph)


first_col_to_rowname <-  function(df){
  rownames(df) <-  df[,1]
  df <-  df[,-1, drop = F]
  return(df)
}

StrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

integrated_for_cor_and_pcor <- function(
    left_ls,right_ls,confounder_ls = NULL, 
    no_pcor_Test_code){ 

  ### 相关分析
  for(i in 1:length(left_ls$output$final_data)){
    if(i == 1) {
      final_ls_for_combination <-  list()
      final_ls_for_combination$Have_partial_result <-  F
      final_df <-  data.frame(stringsAsFactors = F)
    }
    for(j in 1:length(right_ls$output$final_data)){

      if(left_ls$output$final_group[[i]] != right_ls$output$final_group[[j]]){next()}
      cat(i,j,"is being processed.\n")
      left  <-  left_ls$output$final_data[[i]] %>% t()
      right <-  right_ls$output$final_data[[j]] %>% t()
      spearman_df <- spearman_between_two_df_enhanced(left = left,right = right)
      spearman_df <- mutate(spearman_df,
                            pair_1 = pair_1 %>% as.character(),
                            pair_2 = pair_2 %>% as.character() )
      
      right_name <-  names(right_ls$output$final_data)[j] 
      left_name  <-  names(left_ls$output$final_data)[i]  %>% sub("^[^#]+##"," ",.)
      
      ### 核心代码： 将结果保存到列表中，后面再进行可视化操作。
      each_group_info <- left_ls$output$final_group[[i]]
      each_left <- left_name  
      each_right <-  right_name  %>% sub("^[^#]+##"," ",.)
      
      final_ls_for_combination[[each_left]][[each_right]][[each_group_info]][["spearman"]] <- spearman_df
      final_ls_for_combination[[each_left]][[each_right]][[each_group_info]][["right_name"]] <- right_name 
      final_ls_for_combination[[each_left]][[each_right]][[each_group_info]][["left_name"]] <- left_name
      
      each_line <-  data.frame(each_left , each_right, each_group_info,IS_Partial = "spearman",
                               recommend_folder_name = right_name,
                               recommend_file_name = left_name,
                               stringsAsFactors = F)
      final_df <- rbind(final_df, each_line)
      if( is.null(confounder_ls) == F){ ## 不为空的话，执行后面的操作
        
        if(right_ls$note$Test_Code[j] %in%  no_pcor_Test_code| 
           left_ls$note$Test_Code[j] %in%  no_pcor_Test_code){next()}
        
        final_ls_for_combination$Have_partial_result <-  T
        confounder_index <- match(left_ls$output$final_group[i],confounder_ls$output$final_group) ## 通过匹配的方式寻找对应的偏相关变量
        confounding <- confounder_ls$output$final_data[[confounder_index]] %>% t()
        partial_df <- zd_pcor_enhanced(left = left,right = right,confounding = confounding ) ## 偏相关部分
        partial_df <- mutate(partial_df,
                             Method = Method %>% as.character(),
                             pair_1 = pair_1 %>% as.character(),
                             pair_2 = pair_2 %>% as.character() )
        each_line <-  data.frame(each_left , each_right, each_group_info,IS_Partial = "partial",
                                 recommend_folder_name = right_name,
                                 recommend_file_name = left_name, stringsAsFactors = F)
        final_df <- rbind(final_df, each_line)
        final_ls_for_combination[[each_left]][[each_right]][[each_group_info]][["partial"]] <- partial_df
      }else {
        final_ls_for_combination$Have_partial_result <-  F
      }
    }
  }
  
  final_ls_for_combination$final_df_for_each_test <- final_df
  return(final_ls_for_combination)   
}

return_spearman_result_df <-  function(x,y){
  ## 计算相关
  if( sd(x)== 0 | sd(y) == 0){
    p <-  1
    s_Statistic <- 0
    r <-  0
  }else{
    result <-  cor.test( ~ x + y, method = "spearman", exact = T)
    p <-  result$p.value
    s_Statistic <-  result$statistic
    r   <-  result$estimate ### 相关系数 rho （即为R）
  }
  result_df <-  data.frame(p=p, s_Statistic=s_Statistic, r = r)
  return(result_df)
}

spearman_between_two_df_enhanced <-  function(left=left,right=right){

  core_rownames <- intersect(rownames(left),rownames(right)) %>% sort(.)
  if(length(core_rownames) == 0){
    cat("名称匹配出现问题，请重新核查\n")
  }
  ### 提取出交集的变量
  left  <-  left[core_rownames,,drop= F] %>% data.frame(stringsAsFactors = F,check.names = F)
  right <-  right[core_rownames,,drop = F] %>% data.frame(stringsAsFactors = F,check.names = F)
  ### 进行分析
  for ( i in 1:ncol(left)){
    if(i == 1){pooled_df <-  data.frame(stringsAsFactors = F)}
    for(j in 1:ncol(right)){
      x <-  left[,i]
      y <-  right[,j]
      pair_1 <-  colnames(left)[i]
      pair_2 <-  colnames(right)[j]
      #### 进行相关分析
      temp_df <- return_spearman_result_df (x,y)
      temp_df <-  data.frame(pair_1 = pair_1, pair_2 = pair_2, temp_df)
      pooled_df <-  rbind(pooled_df,temp_df)
    }
  }
  ## 获取其他p值
  pooled_df$BH <-  pooled_df$p %>% p.adjust(.,method = "BH")
  pooled_df$holm <-  pooled_df$p %>% p.adjust(.,method = "holm")
  pooled_df$bonferroni <-  pooled_df$p %>% p.adjust(.,method = "bonferroni")
  
  return(pooled_df)
}


get_colors_and_breaks <- function(max = 3,min = -3,breaks = 0.05,
                                  white_color_expend = 2,
                                  blue_raw,red_raw,zero_for_balance = 0,
                                  balance_color  = "white"){
  
  
  result <- list()
  red_part  <- max - zero_for_balance
  blue_part <- zero_for_balance - min
  
  if(white_color_expend == 0){
    white_color <-  balance_color
  }else{
    white_color <- colorRampPalette(c(blue_raw[length(blue_raw)],balance_color,red_raw[1]))(1+white_color_expend *2)
  }
  red_color   <- colorRampPalette(red_raw)(red_part/breaks-white_color_expend)
  blue_color  <- colorRampPalette(blue_raw)(blue_part/breaks-white_color_expend)
  
  my_color  <- c(blue_color,white_color,red_color) ### order is from lower to upper from left to the right.
  my_breaks <- seq(-1*blue_part,red_part, by = breaks)  
  
  result$my_color  <-  my_color
  result$my_breaks <-  my_breaks
  
  return(result)
}

temp_p_trans <-  function ( input , digits= 4 ,
                            breaks = c( -0.00000000001,0.001,0.01,0.05,1.0001) ,
                            labels = c("**","**","*","")
                            , include_right_bound = F  ){  
  
  input <-  input %>% as.character() %>% as.numeric()
  label_vector   <-  cut ( input , breaks , labels , right = include_right_bound  ) %>% str_pad ( . , width = 1, side = "right") 
  neat_number    <-  neat_round ( input , digits =  digits )
  neat_number [ input <0.0001] <- "<0.0001"
  output <-  paste (label_vector ) 
  return ( output )
}

get_pheatmap_enhanced <-  function(final = spearman_result_df,pheatmap_filename = pheatmap_filename,
                                   r_cutoff = 0.3,
                                   p_cutoff = 0.05,
                                   IS_cluster_rows  = T,## 是否行聚类
                                   IS_cluster_cols  = T,## 是否列聚类
                                   simple_taxonname = T, ## 是否简化菌的名字。
                                   taxonname_split = "|",
                                   blue_raw = c("#4292C6","#6BAED6","#9ECAE1","#C6DBEF","#DEEBF7"),
                                   red_raw = c("#F9EFEF","#FEE0D2",'#FCBBA1',"#FC9272","#FB6A4A"),
                                   p_trans_breaks = c( -0.00000000001,0.001,0.01,0.05,1.0001),
                                   p_trans_labels = c("***","**","*",""),
                                   pic_width = 50,
                                   pic_height = 20
                                   
){
 
  r_df <-  final %>% dcast(., pair_1 ~ pair_2, value.var = "r", fun.aggregate = mean) %>% first_col_to_rowname()
  p_df <-  final %>% dcast(., pair_1 ~ pair_2, value.var = "p", fun.aggregate = mean) %>% first_col_to_rowname()

  if(class(r_df) %in% c("numeric","character","factor")){
    return("变量太少 不适合做热图")
  }
  
  IS_r_qualified <- apply(r_df,1,function(each_row){
    IS_row_qualified <- max(abs(each_row), na.rm = T) > r_cutoff
  })
  IS_p_qualified <- apply(p_df,1,function(each_row){
    IS_row_qualified <- min(each_row, na.rm = T) < p_cutoff
  })
  
  r_df <-  r_df[IS_r_qualified & IS_p_qualified,,drop = F]
  p_df <-  p_df[IS_r_qualified & IS_p_qualified,,drop = F]
  if(simple_taxonname == T){ ## 简化菌的名字
    text_for_gsub <-  paste0(".*[",taxonname_split,"](.*)")
    rownames(r_df) <-  rownames(r_df)  %>% gsub(";$","",.) %>%  gsub(text_for_gsub,"\\1",.)
    colnames(r_df) <-  colnames(r_df)  %>% gsub(";$","",.) %>%  gsub(text_for_gsub,"\\1",.)
  } 
  
  if(nrow(p_df) == 0){
    cat("没有符合条件的，做不了\n")
    return("没有符合条件的，做不了")
  }
  ### 对r、p的矩阵进行筛选。
  p_df_label <- apply(p_df,2,function(each_col) temp_p_trans(each_col,
                                                             breaks = p_trans_breaks,  #breaks = c( -0.00000000001,0.001,0.01,0.05,1.0001)
                                                             labels = p_trans_labels#labels = c("**","**","*","")
  )
  )

  if(class(p_df_label)[1] == "character"){ ### 当仅有一行时，数据会被转换为字符向量，再转换回数据框。
    IS_cluster_rows <-  F ### 不再进行聚类。
    p_df_label <- p_df_label %>% data.frame() %>% t()
  }
  
  
  if(nrow(r_df) == 1){ ### 若数据仅一行，不进行聚类
    IS_cluster_rows <-  F
  }
  
  if(ncol(r_df) == 1){ ### 若数据仅一行，不进行聚类
    IS_cluster_cols <-  F
  }
  
  result <- get_colors_and_breaks(max =max( max(r_df), 0.3),
                                  min =min( min(r_df), -0.3),
                                  blue_raw = blue_raw,
                                  red_raw = red_raw,
                                  breaks = 0.002,
                                  white_color_expend = 4)
  
  pic <-  pheatmap(r_df, display_numbers = p_df_label,
                   color  = result$my_color,
                   breaks = result$my_breaks,
                   border_color = "white",
                   number_color = "white",
                   fontsize_number = 12, 
                   family = "TT",
                   angle_col = 45, cellwidth = 15, cellheight = 15,
                   cluster_cols = IS_cluster_cols,
                   cluster_rows = IS_cluster_rows)
  
  paremeter_codes <- paste("p",p_cutoff,"r",r_cutoff)
  
  pheatmap_filename <- pheatmap_filename %>% sub("(.*)(.pdf$|.tiff$|.png$)","\\1",.) %>% paste(paremeter_codes) %>% paste0(., pheatmap_filename %>% sub("(.*)(.pdf$|.tiff$|.png$)","\\2",.)   )
  
  ggsave(pheatmap_filename, plot = pic , width =pic_width , height = pic_height, units = c("cm"),scale = 1, limitsize = F)
}


getwd()
get_or_set_dir()

index <- 1
micro<- "Microbiota.xlsx"
meta <- "Metabolite.xlsx"


meta_ls <- Get_ls_from_exl_5.0(   ## 代谢数据
  file_name = meta,
  IS_microbiota = F,
  remove_letter__ = F)


micro_ls <- Get_ls_from_exl_5.0( ## 菌群数据
  file_name = micro,
  IS_microbiota = T,
  sep_code =  ";",
  get_level_by = "sep_code",
  remove_letter__ = T,
  match_mode = "tail_match",   ## 匹配方式tail_match,截断匹配
  adj_index =  1
)


{final_ls_for_next <- integrated_for_cor_and_pcor(
  left_ls = micro_ls,  ##  实际画在右边
  right_ls = meta_ls   ##  实际画在左边
) 
  
  Note_for_net <- c("Microbiota","Metabolite") ## 分别对应left 和 right
  
  for(i in 1:nrow(final_ls_for_next$final_df_for_each_test)){
    if(i == 1){
      pooled_df <- data.frame(stringsAsFactors = F, check.names = F)
      index_df <- final_ls_for_next$final_df_for_each_test
      
    }
    
    index <- index_df[i,]
    final_df <- final_ls_for_next[[index$each_left]][[index$each_right]][[index$each_group_info]][[index$IS_Partial]]
    temp_df <- data.frame(group = index$each_group_info,
                          bottom_in_heatmap  = index$each_right,
                          right_in_heatmap  = index$each_left,
                          Method = index$IS_Partial,
                          final_df, check.names = F, stringsAsFactors = F) 

    pooled_df <-  rbind(pooled_df, temp_df)
    final_ls_for_next$final_df_for_each
    
    sub1 <- micro %>% sub("\\.xlsx","",.)
    sub2 <- meta %>% sub("\\.xlsx", "", .)
    path <-  paste0(output_path, index$recommend_folder_name, sub1,"#", sub2, "/")
    dir.create(path,recursive = T, showWarnings = F)
    file_trunk <- paste0(path,index$IS_Partial, index$recommend_file_name)
    file_name <- paste0(file_trunk,'_', sub1,'_',sub2, ".csv")
    
    heatmap_name <- paste0(file_trunk," heatmap.pdf")

    
    pic_height <- length(unique(final_df$pair_1))*0.6+20
    pic_width <- length(unique(final_df$pair_2))*0.6+20
    write.csv(final_df, file_name)  ### 写入文件
    get_pheatmap_enhanced(final = final_df,
                          pheatmap_filename = heatmap_name,
                          pic_width = pic_width,
                          pic_height =  pic_height,
                          simple_taxonname = T,
                          r_cutoff = 0,
                          p_cutoff = 1,
                          taxonname_split = ";")
    
    get_pheatmap_enhanced(final = final_df,
                          pheatmap_filename = heatmap_name,
                          pic_width = pic_width,
                          pic_height =  pic_height,
                          simple_taxonname = T,
                          r_cutoff = 0.3,
                          p_cutoff = 0.05,
                          taxonname_split = ";")
    
    get_pheatmap_enhanced(final = final_df,
                          pheatmap_filename = heatmap_name,
                          pic_width = pic_width,
                          pic_height =  pic_height,
                          simple_taxonname = T,
                          r_cutoff = 0.5,
                          p_cutoff = 0.05,
                          taxonname_split = ";")
    
    get_pheatmap_enhanced(final = final_df,
                          pheatmap_filename = heatmap_name,
                          pic_width = pic_width,
                          pic_height =  pic_height,
                          simple_taxonname = T,
                          r_cutoff = 0.7,
                          p_cutoff = 0.05,
                          taxonname_split = ";")
  }
}
