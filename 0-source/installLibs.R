installLibs <- function(libs,
                        type = 'CRAN'){
  
  for (i in libs){
    
    if(type != 'github'){
    cat('> package ', i, sep = '')
    } else {
      cat('> package ', names(libs)[which(libs==i)], sep = '')
    }
    
    installed <- rownames(installed.packages())
    
    if(type == 'CRAN'){
      
      if(!i %in% installed){
        cat(': installing ...')
        install.packages(i, character.only = T,quiet = T)
        cat('done\n')
      } else {
        cat(': already installed\n')
      }
    }
    
    if(type == 'bioc'){
      
      installed <- rownames(installed.packages())
      if (!"BiocManager" %in% installed){
        install.packages("BiocManager")
        }
      
      if(!i %in% installed){
        cat(': installing ...')
        BiocManager::install(i, character.only = T,ask = F)
        cat('done\n')
      } else {
        cat(': already installed\n')
      }
      
    }
    
    
    if(type == 'github'){
      
      if(!'devtools' %in% installed){install.packages("devtools", quiet = T)}
      
      if(!names(libs)[which(libs==i)] %in% installed){
        cat(': installing ...')
        devtools::install_github(i)
        cat('done\n')
      } else {
        cat(': already installed\n')
      }
      
    }
    
  }
  
  cat("All ", type, " libraries installed.\n\n", sep = "")
  
}
