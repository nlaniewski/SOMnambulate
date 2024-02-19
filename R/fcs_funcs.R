#' @title Get parameters ('$P') from .fcs file headers
#' @description
#' $P keywords in the TEXT section of a .fcs file denote cytometer-specific parameter names/values. There is useful information that can be parsed from these parameters -- particularly the $PN and $PS parameter-value pairs. This function allows parsing of these parameters without the overhead of reading the actual .fcs file to get the values.
#'
#'
#' @param fcs.file.paths Character string; path(s) usually returned from \code{list.files(...,full.names=T,pattern=".fcs")}.
#' @param return.dt Logical. By default, \code{FALSE}; if \code{TRUE}, will return a data.table of parameters per .fcs file
#'
#' @return a list of parameters per .fcs file; if \code{return.dt=T}, a data.table of parameters
#' @examples
#' extdata<-system.file("extdata",package="SOMnambulate")
#' fcs.files <- list.files(extdata,full.names=TRUE,pattern=".fcs")
#' get.fcs.parameters(fcs.files[1])
#' get.fcs.parameters(fcs.files[1],return.dt=TRUE)[]
#'
#' @export
get.fcs.parameters<-function(fcs.file.paths,return.dt=F){
  fcs.parameters.list<-sapply(flowCore::read.FCSheader(fcs.file.paths),function(h){
    p.max<-as.numeric(h['$PAR'])
    parm.letters<-unique(sub("\\$P[0-9]+","",grep('\\$P[0-9]',names(h),value = T)))
    p.list<-sapply(parm.letters,function(p){
      pars<-h[grep(paste0("\\$P[0-9]+",p),names(h))]
    },simplify = F)
    for(i in names(which(sapply(p.list,length)<p.max))){
      fill<-setdiff(seq(p.max),as.numeric(stringr::str_extract(names(p.list[[i]]),"[0-9]+")))
      fill.names<-stats::setNames(nm=paste0("$P",fill,i),rep(NA,length(fill)))
      p.list[[i]]<-c(p.list[[i]],fill.names)
    }
    p.list<-lapply(p.list,function(p) p<-p[stringr::str_order(names(p),numeric = T)])
    return(p.list)
  },simplify = F)
  if(return.dt){
    fcs.parameters.list<-sapply(fcs.parameters.list,function(p) data.table::setDT(p),simplify = F)
  }else{
    return(fcs.parameters.list)
  }
}
#' @title Get keyword metadata (!'$P') from .fcs files
#' @description
#' Non-parameter ($P) keywords in the TEXT section of a .fcs file denote cytometer/sample-specific meta-data. This function allows parsing of that meta-data without the overhead of reading the actual .fcs file to get the values.  Depending on the nature of the .fcs file, there will be some mix of FCS-defined ('$') and non-standard keywords.
#'
#'
#' @param fcs.file.paths Character string; path(s) usually returned from \code{list.files(...,full.names=T,pattern=".fcs")}.
#' @param return.dt Logical. By default, \code{FALSE}; if \code{TRUE}, will return a data.table of keyword/values per .fcs file
#' @param pattern Pattern used to split keyword character strings; specifically for use with barcoded/pooled .fcs files that contain 'collapsed' metadata.
#' @param pattern.unique Logical. By default, \code{TRUE}; returns unique keyword data.tables.
#'
#' @return a list of non-parameter ('$P') keywords per .fcs file; if \code{return.dt=T}, a data.table of keywords/values
#' @examples
#' extdata<-system.file("extdata",package="SOMnambulate")
#' fcs.files <- list.files(extdata,full.names=TRUE,pattern=".fcs")
#' get.fcs.keywords.metadata(fcs.files[1])
#'
#' #note column classes in the data.table; the TEXT section is encoded as character
#' get.fcs.keywords.metadata(fcs.files[1],return.dt=TRUE)[]
#'
#' #list; row-bound data.table; a few example conversions
#' mdats <- get.fcs.keywords.metadata(fcs.files,return.dt=TRUE)#a list of data.tables
#' mdat<-data.table::rbindlist(mdats)#a single, row-bound data.table
#' mdat[,class(`$DATE`)]#"character" class
#' mdat[,`$DATE`:= data.table::as.IDate(`$DATE`,format="%d-%b-%Y")]
#' mdat[,class(`$DATE`)]#"Date" class
#'
#' @export
get.fcs.keywords.metadata <- function(fcs.file.paths,return.dt=F,pattern=NULL,pattern.unique=T){
  fcs.keywords.list <- sapply(flowCore::read.FCSheader(fcs.file.paths),function(h){
    #drop-terms: '$P[0-9]+' (parameters); 'P[0-9]+DISPLAY' (Cytek specific?); Spill(over)
    h<-h[-grep('\\$PAR|\\$P[0-9]+|P[0-9]+DISPLAY|spill',names(h),ignore.case = T)]
    return(as.list(h))
  },simplify = F)
  if(return.dt&is.null(pattern)){
    fcs.keywords.dt <- sapply(fcs.keywords.list, function(kw) data.table::setDT(kw),simplify = F)
  }else if(return.dt&!is.null(pattern)){
    fcs.keywords.dt.split <- lapply(fcs.keywords.list,function(kws){
      lapply(sort(unique(stringr::str_count(kws,pattern))),function(s){
        data.table::as.data.table(sapply(kws[which(stringr::str_count(kws,pattern)==s)],strsplit,pattern))
      })
    })
    if(pattern.unique){
      fcs.keywords.dt.split<-sapply(seq(max(sapply(fcs.keywords.dt.split,length))),function(n){
        fcs.keywords.dt.split<-sapply(fcs.keywords.dt.split,'[',n)
        if(length(unique(fcs.keywords.dt.split))==1) fcs.keywords.dt.split<-fcs.keywords.dt.split[[1]]
        return(fcs.keywords.dt.split)
      })
    }else{
      return(fcs.keywords.dt.split)
    }
  }else{
    return(fcs.keywords.list)
  }
}
#' @title Get a `data.table` of .fcs file paths and related metadata.
#' @description
#' Essentially a wrapper around \link{get.fcs.keywords.metadata} with some convenience formatting and class conversions. This structured `data.table` can be considered the starting point of an analysis workflow. Naming conventions can be enforced by modifying any sample/subject-specific columns; it's not uncommon for typos to occur during the course of a large study and they can be corrected after the keyword metadata has been parsed into this `data.table`.
#'
#'
#' @param fcs.file.paths Character string; path(s) usually returned from `list.files(...,full.names=T,pattern=".fcs")`.
#' @param factor.cols Character string; column names to be converted to factor.
#'
#' @return a `data.table` of full length .fcs file paths and related metadata (parsed from the TEXT header).
#' @examples
#' extdata<-system.file("extdata",package="SOMnambulate")
#' fcs.files <- list.files(extdata,full.names=TRUE,pattern=".fcs")
#' fcs.files.dt<-get.fcs.file.dt(fcs.files)
#' fcs.files.dt[]
#'
#' #enforce/establish naming convention for the construction of a sample.id and associated factors
#' new.cols<-c('study.name','batch.seq','batch.date')
#' fcs.files.dt[,(new.cols) :=
#' data.table::tstrsplit(basename(dirname(FILENAME)),"_",type.convert=list(as.factor=1:3))]
#' fcs.files.dt[,stim.condition := factor(stringr::str_extract(FIL,"SEB|UNSTIM"))]
#' fcs.files.dt[,aliquot.seq := factor(seq(.N)),by=.(batch.date,stim.condition)]
#' fcs.files.dt[,batch := factor(paste(study.name,batch.seq,batch.date,sep="_"))]
#' fcs.files.dt[,sample.id := paste(batch,stim.condition,aliquot.seq,sep="_")]
#' fcs.files.dt[]
#'
#' @export
get.fcs.file.dt<-function(fcs.file.paths,factor.cols=NULL){
  f.path<-source.name<-file.size.MB<-NULL
  dt<-data.table::data.table(f.path=fcs.file.paths)
  dt[,source.name:=basename(f.path)]
  dt[,file.size.MB:=signif(file.size(f.path)/1024^2,3)]
  ##
  dts<-get.fcs.keywords.metadata(fcs.file.paths,return.dt = T)
  dt<-cbind(dt,data.table::rbindlist(dts,fill=T))
  #synatically valid names; drop keyword identifier '$'
  names(dt)<-sub("\\$","",names(dt))
  #drop conserved/standard keywords that are non-informative to the user
  drop.terms<-c(
    paste0(c('BEGIN','END'),c('ANALYSIS')),
    paste0(c('BEGIN','END'),c('DATA')),
    paste0(c('BEGINS','ENDS'),c('TEXT')),
    "MODE",
    'BYTEORD',
    'DATATYPE',
    'NEXTDATA',
    'TIMESTEP',
    'APPLY COMPENSATION',
    'USERSETTINGNAME',
    'CHARSET',
    'TargetGemStoneMethod',#mass cytometry-specific; Helios (CYTOF-7.0.5189)
    'FluidigmBarcoded',#mass cytometry-specific; Helios (CYTOF-7.0.5189)
    'FluidigmMassTagsChecked'#mass cytometry-specific; Helios (CYTOF-7.0.5189)
  )
  dt<-dt[,!drop.terms[drop.terms %in% names(dt)],with=F]
  #do a few conversions
  for(j in c('DATE')){data.table::set(dt,j=j,value=data.table::as.IDate(dt[[j]],format="%d-%b-%Y"))}
  for(j in c('BTIM','ETIM')){data.table::set(dt,j=j,value=data.table::as.ITime(dt[[j]]))}
  for(j in grep("TOT|DELAY|ASF|VOL$",names(dt),value = T)){data.table::set(dt,j=j,value=as.numeric(dt[[j]]))}
  #conserved keywords -- convert these to factor
  keywords.conserved <- c('CYT','CYTSN','FCSversion')
  for(j in keywords.conserved){data.table::set(dt,j=j,value=factor(dt[[j]]))}
  #non-standard, mass cytometry specific -- DVSSCIENCES-FLUIDIGM-CYTOF-7.0.5189; logical
  keywords.fluidigm <- c('FluidigmNormalized')
  for(j in keywords.fluidigm){data.table::set(dt,j=j,value=as.logical(dt[[j]]))}
  #defined through argument
  if(!is.null(factor.cols)){
    for(j in factor.cols){if(j %in% names(dt)) data.table::set(dt,j=j,value=factor(dt[[j]]))}
  }
  ##
  return(dt)
}
#' @title Get a `channel_alias` data.table from .fcs file headers
#' @description
#' The resultant `data.table` is to be used with the `channel_alias` argument of \link[flowCore]{read.FCS}. The alias will be used to set column names; make sure they are unique and 'readable'. For example, the included mass cytometry data files used in this example are originally named ($PS) as 'MassMetal_Marker'; a more 'readable' form would be 'Marker_MassMetal' or even just 'Marker' (see example). The intent here is to use unified names that follow an established convention. There is relatively low 'overhead' reading .fcs headers/TEXT, so having a unified 'channel.alias' (panel) before reading in often large amounts of .fcs data will streamline the workflow.
#'
#' @param fcs.file.paths Character string; path(s) usually returned from `list.files(...,full.names=T,pattern=".fcs")`.
#' @param drop.pattern 'pattern' argument for \link[base]{grep}: "character string containing a \link[base:regex]{regular expression}"; will drop pattern matches from the alias column of `channel_alias`.
#' @param name.sub Named character vector for use in resolving name conflicts/discrepancies; the vector element(s) should equal a pattern and the name(s) a replacement string.
#' @param order.alias Logical; if `TRUE`, the `channel_alias` 'alias' column will be sensibly ordered.
#'
#' @return returns a `data.table` containing a 'channels' and 'alias' column; returns a list of `data.table`s if not unique (issues a warning).
#' @examples
#' extdata<-system.file("extdata",package="SOMnambulate")
#' fcs.files <- list.files(extdata,full.names=TRUE,pattern=".fcs")
#'
#' #a data.table containing 'channels' ($PN) and 'alias' ($PS)
#' get.fcs.channel.alias(fcs.files)
#'
#' #extraneous channels can be dropped
#' drop.metals<-c("89Y",paste0(c(102,104:106,108,110),"Pd"),"120Sn","127I",
#' paste0(c(138,139),"La"),"176Lu","190Os",paste0(c(191,193),"Ir"),
#' paste0(c(195,198),"Pt"),"208Pb")
#' drop.metals<-paste0(drop.metals,'$',collapse="|")
#' drop.pattern<-paste(drop.metals,'background','noise',sep="|")
#'
#' ca<-get.fcs.channel.alias(fcs.files,drop.pattern=drop.pattern,order.alias=TRUE)
#' ca[]
#'
#' #a few more 'name-fixes'
#' ca[grep("beads",alias,ignore.case = TRUE),alias:=sub("Norm_beads","beads",alias)]
#' ca[grep("viability",alias,ignore.case = TRUE),alias:="194Pt_viability"]#drop '_cisplatin'
#' ca[]
#' #flip metal and marker to make more 'readable'
#' ca[stringr::str_detect(alias,"[A-Z]{1}[a-z]{1}_"),alias := sub('(\\w+)_(\\w+)', '\\2_\\1', alias)]
#' ca[]
#'
#' #for the sake of this example, assume a .fcs parameter was misnamed
#' fcs.tmp<-flowCore::read.FCS(fcs.files[1],transformation=FALSE,truncate.max.range=FALSE)
#' flowCore::markernames(fcs.tmp)[['Dy164Di']]<-'164Dy_IFNy'
#' fcs.tmp.path<-tempfile(fileext = ".fcs")
#' flowCore::write.FCS(fcs.tmp,fcs.tmp.path)
#'
#' #a warning will issue; two unique data.tables will be returned;
#' get.fcs.channel.alias(fcs.file.paths=c(fcs.files[2],fcs.tmp.path),
#' drop.pattern=drop.pattern)
#'
#' #using name.sub argument to resolve the warning and return a single data.table
#' get.fcs.channel.alias(fcs.file.paths=c(fcs.files[2],fcs.tmp.path),
#' drop.pattern=drop.pattern, name.sub=c("164Dy_IFNg"="164Dy_IFNy"))
#'
#' @export
get.fcs.channel.alias<-function(fcs.file.paths,drop.pattern=NULL,name.sub=NULL,order.alias=F){
  channel_alias.list<- unique(lapply(fcs.file.paths,function(f.path){
    ca<-SOMnambulate::get.fcs.parameters(f.path,return.dt = T)[[1]]
    ca[is.na(S), S := N]
    if(ca[,data.table::uniqueN(S)]!=ca[,.N]){
      stop("Non-unique alias name")
    }else{
      data.table::setnames(ca,c('N','S'),c('channels','alias'))
      if(!is.null(drop.pattern)){
        ca<-ca[grep(drop.pattern,alias,ignore.case=TRUE,invert=TRUE),.(channels,alias)]
      }else{
        ca<-ca[,.(channels,alias)]
      }
    }
    if(!is.null(name.sub)){
      for(n in names(name.sub)){
        data.table::set(ca,i=grep(name.sub[[n]],ca$alias),j='alias',value=n)
      }
    }
    if(order.alias){
      if(any(grepl("event_length",ca$channels,ignore.case = T))){
        alias.order<-c('Time',sort(grep('event|center|offset|width|residual',ca$alias,value = T,ignore.case = T)))#mass cytometry
      }else if(any(grepl("FSC|SSC",ca$channels,ignore.case = T))){
        alias.order<-c('Time',sort(grep('FSC|SSC',ca$alias,value = T,ignore.case = T)))#conventional/spectral cytometry
      }
      alias.order<-c(alias.order,ca[!alias %in% alias.order,stringr::str_sort(alias,numeric = T)])
      data.table::setorder(ca[, 'ord' := match(alias,alias.order)],'ord')[,'ord' := NULL]
    }
    return(ca)
  }))
  ##
  if(length(channel_alias.list)==1){
    return(channel_alias.list[[1]])
  }else{
    warning(paste(
      'channel and/or alias conflict; resolve name discrepancy:',
      paste0(setdiff(
        Reduce(union,lapply(channel_alias.list,'[[','alias')),
        Reduce(intersect,lapply(channel_alias.list,'[[','alias'))),collapse = " : ")
    ),call. = F)
    return(channel_alias.list)
  }
}
#' @title Convert a \code{data.table} of mass cytometry .fcs data into a new .fcs file
#' @description
#' After reading in (usually large amounts of) .fcs data and converting to a \code{data.table}, the process can be reversed and a new .fcs file created; new keywords are written to the header section and parameters updated using \code{flowCore} functions.
#'
#'
#' @param dt.data A \code{data.table} of .fcs data, as returned from \code{data.table::as.data.table(flowCore::read.FCS(...))}.
#' @param reverse.asinh.cofactor If defined, will reverse the \code{asinh} transformation - using the same cofactor value - typically applied to mass cytometry .fcs data
#' @param keywords.to.add A list of keywords to append to the header section of the .fcs file; if \code{NULL}, returns a 'basic' \code{flowCore::flowFrame}.
#' @param fil If defined, will update the '$FIL' keyword; when writing the new .fcs file, '$FIL' is used to define the \code{flowCore::write.FCS(filename)} argument.
#' @param out.path A \code{file.path}; will write the new .fcs to the defined directory; if \code{NULL}, returns an 'updated' \code{flowCore::flowFrame}.
#'
#' @return Depending on defined arguments, returns a \code{flowCore::flowFrame} with/without updated keywords or writes the .FCS to a defined directory.
#' @export
#'
fcs.from.dt.masscyto<-function(dt.data,reverse.asinh.cofactor=NULL,keywords.to.add=NULL,fil=NULL,out.path=NULL){
  maxRange<-NULL
  if(!is.null(reverse.asinh.cofactor)){
    if(!is.numeric(reverse.asinh.cofactor)) stop("Non-numeric value defined for 'reverse.asinh.cofactor'")
    cols.metal <- grep("[0-9]{3}[A-Z]{1}[a-z]{1}",names(dt.data),value = T)
    dt.data[ , (cols.metal) := lapply(.SD,function(x){sinh(x)*reverse.asinh.cofactor}), .SDcols = cols.metal]
  }
  parms.adf<-data.table::data.table(name = names(dt.data),desc = NA,dt.data[,.(minRange=0,maxRange=unlist(lapply(.SD,function(x) ceiling(max(x)))))])
  parms.adf[,range:=(maxRange+1)]
  ##
  parms.list<-as.list(
    stats::setNames(c(rep(c("32","0,0"),each=nrow(parms.adf)),parms.adf$range),
                    nm=do.call(paste0, expand.grid(paste0('$P',seq(nrow(parms.adf))), c("B","E","R"))))
  )
  ##
  if(!is.null(keywords.to.add)){
    fcs.from.dt<-methods::new("flowFrame",exprs=as.matrix(dt.data),parameters=Biobase::AnnotatedDataFrame(parms.adf),description = parms.list)
    flowCore::keyword(fcs.from.dt)<-c(flowCore::keyword(fcs.from.dt),keywords.to.add)
    flowCore::keyword(fcs.from.dt)[['modified.by']]<-"R_fcs.from.dt.masscyto"
    if(!is.null(fil)){
      flowCore::keyword(fcs.from.dt)[['$FIL']] <- fil
    }
    if(is.null(out.path)){
      return(fcs.from.dt)
    }else{
      if(!dir.exists(out.path)) dir.create(out.path,recursive = T)
      flowCore::write.FCS(fcs.from.dt,file.path(out.path,flowCore::keyword(fcs.from.dt)[['$FIL']]))
    }
  }else{
    return(methods::new("flowFrame",exprs=as.matrix(dt.data),parameters=Biobase::AnnotatedDataFrame(parms.adf),description = parms.list))
  }
}
#' @title reads a .fcs file and converts it to `data.table`
#' @description
#' After reading a .fcs file, the expression matrix is converted to a `data.table`; in addition, character and factor columns (meta-data) are appended to the numeric expression values. Once in `data.table` form, operations such as transformation and scaling can be done 'in place' and 'by' groupings streamline analyses. The example data size is fractionally small compared to even a single typical .fcs file; with sufficient RAM resources, a `data.table` of .fcs data amounting to tens of GBs (many high-parameter (spectral/mass) .fcs files, row-bound) can be held in memory.
#'
#'
#' @param fcs.file.dt a `data.table` as returned from \link{get.fcs.file.dt}; `fcs.file.dt` can be subset to include only columns-of-interest (a minimally-representative subset would include only 'f.path' and 'sample.id').
#' @param channel_alias as returned from \link{get.fcs.channel.alias}
#' @param use.alias.pattern Logical; default `FALSE`. If `TRUE` and `channel_alias` is defined, the 'alias' column will be used as a pattern to define the \link[flowCore]{read.FCS} `column.pattern` argument.
#' @param alias.order Logical. If `TRUE`, the `data.table` columns will be ordered to match that of the `channel_alias` 'alias' column.
#' @param cofactors A named numeric vector; named columns will be \link[base]{asinh} transformed with the supplied cofactor (numeric).
#'
#' @return a `data.table` of raw, un-transformed numeric expression values with character/factor identifier columns; if `cofactors` is defined, the raw expression values will be \link[base]{asinh} transformed.
#' @examples
#' \dontrun{
#' #from the 'get.fcs.file.dt' example:
#' fcs.files.dt_ECHO
#'
#' #from the 'get.fcs.channel.alias' example:
#' ca_ECHO
#'
#' #retains original 'detector' names;
#' #all columns from 'fcs.file.dt'
#' names(SOMnambulate:::fcs.to.dt(fcs.files.dt_ECHO[1],channel_alias=NULL))
#'
#' #replaces the original names with the 'alias' from 'channel.alias';
#' #all columns from 'fcs.file.dt'
#' names(fcs.to.dt(fcs.files.dt_ECHO[1],channel_alias=ca_ECHO))
#'
#' #replaces the original names with the 'alias' from 'channel.alias'; drops via pattern
#' #all columns from 'fcs.file.dt'
#' names(fcs.to.dt(fcs.files.dt_ECHO[1],channel_alias=ca_ECHO,use.alias.pattern=TRUE))
#'
#' #replaces the original names with the 'alias' from 'channel.alias'; drops via pattern
#' #reorders columns using 'alias.order'
#' #a subset of columns from 'fcs.file.dt'
#' names(fcs.to.dt(fcs.files.dt_ECHO[1,.(f.path,sample.id,batch,stim.condition,aliquot.seq)],
#' channel_alias=ca_ECHO,use.alias.pattern=TRUE,alias.order=TRUE))
#'
#' dt<-fcs.to.dt(fcs.files.dt_ECHO[1,.(f.path,sample.id,batch,stim.condition,aliquot.seq)],
#' channel_alias=ca_ECHO,use.alias.pattern=TRUE,alias.order=TRUE)
#' dt[]
#'
#' #as a list of individual .fcs files
#' dts<-sapply(
#' split(fcs.files.dt_ECHO[,.(f.path,sample.id,batch,stim.condition,aliquot.seq)],by='f.path'),
#' fcs.to.dt,channel_alias=ca_ECHO,use.alias.pattern=TRUE,alias.order=TRUE,simplify=FALSE)
#' dts[1]
#' }
#'
#' @export
fcs.to.dt<-function(fcs.file.dt,channel_alias=NULL,use.alias.pattern=FALSE,alias.order=FALSE,cofactors=NULL){
  #read .fcs file ('.path')
  #if defined, rename columns using 'channel_alias' as returned from 'get.fcs.channel.alias'
  fcs.tmp<-flowCore::read.FCS(fcs.file.dt[['f.path']],transformation = F,truncate_max_range = F,
                              channel_alias = if(!is.null(channel_alias)) channel_alias,
                              column.pattern = if(!is.null(channel_alias) & use.alias.pattern) paste0(channel_alias$alias,collapse = "|")
  )
  #convert the expression matrix (raw, un-transformed data values) into a data.table
  dt<-data.table::setDT(as.data.frame(fcs.tmp@exprs))
  #if defined, transform named columns with supplied cofactors
  if(!is.null(cofactors)){
    if(!all(names(cofactors) %in% names(dt))){
      not.found<-names(cofactors)[!names(cofactors) %in% names(dt)]
      stop(paste("The following supplied cofactors are not named columns in 'fcs.file.dt':",
                 paste0(not.found,collapse = " ; "),
                 sep = "\n"
                 )
      )
    }else{
      for(j in names(cofactors)){
        data.table::set(dt,j=j,value = asinh(dt[[j]]/cofactors[[j]]))
      }
    }
  }
  #reorder columns according to an ordered alias
  if(alias.order){
    data.table::setcolorder(dt,channel_alias$alias)
  }
  #column-bind identifier columns
  if(length((fcs.file.dt[,!'f.path']))>0){
    dt<-cbind(dt,fcs.file.dt[,!'f.path'])
  }
  #
}
#' @title a parallelized version of \link{fcs.to.dt}
#' @description
#' Essentially a wrapper around \link[parallel]{makeCluster} and those associated functions.  There is 'computational' overhead in initiating the compute clusters but for many .fcs files, this parallelized version of \link{fcs.to.dt} should speed up the workflow.  There is risk of saturating all available RAM if the data set is large enough;  it is recommended to monitor system resources (Task Manager/Activity Monitor).
#'
#' @param fcs.file.dt a `data.table` as returned from \link{get.fcs.file.dt}; `fcs.file.dt` can be subset to include only columns-of-interest (a minimally-representative subset would include only 'f.path' and 'sample.id').
#' @param channel_alias as returned from \link{get.fcs.channel.alias}
#' @param use.alias.pattern Logical; default `FALSE`. If `TRUE` and `channel_alias` is defined, the 'alias' column will be used as a pattern to define the \link[flowCore]{read.FCS} `column.pattern` argument.
#' @param alias.order Logical. If `TRUE`, the `data.table` columns will be ordered to match that of the `channel_alias` 'alias' column.
#' @param cofactors A named numeric vector; named columns will be \link[base]{asinh} transformed with the supplied cofactor (numeric).
#'
#' @return a `data.table` of raw, un-transformed numeric expression values (row-bound) with character/factor identifier columns; if `cofactors` is defined, the raw expression values will be \link[base]{asinh} transformed.
#' @examples
#' \dontrun{
#' #because of the small file-size of these example .fcs files,
#' #it's faster to read them individually with 'fcs.to.dt';
#' #for larger/many files, the overhead of initiating the compute clusters is negligible
#' #compared to the speed savings of parallel reading
#'
#' #from the 'get.fcs.file.dt' example:
#' fcs.files.dt_ECHO
#'
#' #from the 'get.fcs.channel.alias' example:
#' ca_ECHO
#'
#' dt<-fcs.to.dt.parallel(
#' fcs.file.dt=fcs.files.dt_ECHO[,.(f.path,sample.id,batch,stim.condition,aliquot.seq)],
#' channel_alias=ca_ECHO,use.alias.pattern=TRUE,alias.order=TRUE)
#' dt[]
#' }
#'
#' @export
fcs.to.dt.parallel<-function(fcs.file.dt,channel_alias=NULL,use.alias.pattern=FALSE,alias.order=FALSE,cofactors=NULL){
  ##for posterity, leaving the commented code chunk; the function arguments at one point had be exported;
  ##now working as intended in the uncommented code chunk, which used to not work...
  # n.cores<-parallel::detectCores()-1
  # n.paths<-fcs.file.dt[,.N]
  # n<-ifelse(n.paths>n.cores,n.cores,n.paths)
  # cl<-parallel::makeCluster(n)
  # if(!is.null(channel_alias)){
  #   parallel::clusterExport(cl,'channel_alias',envir = environment())
  # }
  # parallel::clusterExport(cl,'use.alias.pattern',envir = environment())
  # parallel::clusterExport(cl,'alias.order',envir = environment())
  # if(!is.null(cofactors)){
  #   parallel::clusterExport(cl,'cofactors',envir = environment())
  # }
  # ##
  # dt<-data.table::rbindlist(parallel::parLapply(cl,split(fcs.file.dt,by='f.path'),fcs.to.dt,
  #                                               channel_alias = if(!is.null(channel_alias)) channel_alias,
  #                                               use.alias.pattern=use.alias.pattern,
  #                                               alias.order=alias.order,
  #                                               cofactors = if(!is.null(cofactors)) cofactors))
  # ##
  # on.exit({parallel::stopCluster(cl);invisible(gc())})
  ##...
  n.cores<-parallel::detectCores()-1
  n.paths<-fcs.file.dt[,.N]
  n<-ifelse(n.paths>n.cores,n.cores,n.paths)
  message(paste("Parallel reading of",n.paths,".fcs files using",n,"cores."))
  cl<-parallel::makeCluster(n)
  on.exit({parallel::stopCluster(cl);invisible(gc())})
  data.table::rbindlist(parallel::parLapply(cl,split(fcs.file.dt,by='f.path'),
                                            fcs.to.dt,
                                            channel_alias=channel_alias,
                                            use.alias.pattern=use.alias.pattern,
                                            alias.order=alias.order,
                                            cofactors=cofactors))
}
#' @title Generate a `data.table` of .fcs parameters
#'
#' @param dt `data.table` returned from `fcs.to.dt`; some mix of .fcs expression values and character/factor values.
#' @param type.convert Character vector; will attempt to convert each string element (named column in `dt`) to numeric
#' @param name.fix a `data.table` with two columns: 'name' and 'name.fix'; if 'name' matches, the value will be replaced with 'name.fix'.
#'
#' @return a `data.table` of .fcs parameters; used to define `flowCore::flowFrame(...,parameters = Biobase::AnnotatedDataFrame(...))`.
#'
#'
#'
parms.adf.dt<-function(dt,type.convert=NULL,name.fix=NULL){
  #for R CMD check; data.table vars
  #name<-desc<-range<-minRange<-maxRange<-value<-NULL
  #
  col.classes<-sapply(dt,class)
  if(is.null(type.convert)&any(unique(col.classes)!='numeric')){
    message(
      paste("Using numeric columns only; use the 'type.convert' argument if the following are to be included:",
            paste("Non-numeric cols:",paste0(names(col.classes[col.classes!="numeric"]),collapse=" ; ")),
            sep = "\n")
    )
  }else if(!is.null(type.convert)){
    type.test<-suppressWarnings(dt[,.SD,.SDcols = type.convert][,lapply(.SD,function(x){any(is.na(as.numeric(x)))})])
    if(any(type.test)){
      stop("Column(s) as defined by 'type.convert' could not be converted to numeric; NAs introduced by coercion")
    }else{
      col.classes[names(col.classes) %in% type.convert]<-'numeric'
    }
  }
  numeric.cols<-names(col.classes)[col.classes=="numeric"]
  #
  parms.adf<-data.table::melt(dt[,.SD,.SDcols = numeric.cols][,lapply(.SD,as.numeric)],
                              measure.vars=numeric.cols,variable.name='name')[
                                ,stats::setNames(as.list(range(value)),nm=c('minRange','maxRange')),by=name]
  parms.adf[,minRange := floor(minRange)];parms.adf[,maxRange := ceiling(maxRange)]
  parms.adf[,range:=(maxRange)-(minRange)]
  if(!is.null(name.fix)){
    parms.adf<-merge(parms.adf,name.fix,all.x = T,sort=F)
    parms.adf[is.na(name.fix),name.fix:=name]
    parms.adf[,name:=name.fix];parms.adf[,name.fix:=NULL]
  }
  #
  parms.adf[,desc:=NA]
  #
  data.table::setcolorder(parms.adf,c('name','desc','range','minRange','maxRange'))
  #
  return(parms.adf[])
}
#' @title Generate a list of .fcs meta-data
#'
#' @param parms.adf `data.table` returned from `parms.adf.dt`
#' @param cyto.method A single argument for now; 'spectral'; depending on cytometry platform -- 'spectral','conventional','mass' -- list types might change.
#'
#' @return a list of .fcs meta-data; used to define `flowCore::flowFrame(...,description = ...)`.
#'
#'
#'
parms.list.from.adf<-function(parms.adf,cyto.method='spectral'){
  #for R CMD check; data.table vars
  name<-desc<-NULL
  #
  parms.list <- as.list(
    c(stats::setNames(
      c(rep(c("32", "0,0"),each = parms.adf[,.N]), parms.adf$range),
      nm = do.call(paste0,expand.grid(paste0("$P", seq(nrow(parms.adf))),c('B','E','R')))))
  )
  parms.list<-append(parms.list,stats::setNames(parms.adf[,name],nm=paste0("$P",seq(parms.adf[,.N]),'N')))
  if(parms.adf[!is.na(desc),.N]>0){
    parms.list<-append(parms.list,stats::setNames(parms.adf[!is.na(desc),desc],nm=paste0("$P",parms.adf[,.I[!is.na(desc)]],'S')))
  }
  #
  if(cyto.method=='spectral'){
    types<-parms.adf[,ifelse(grepl('Time',name),'Time',
                             ifelse(grepl('FSC',name),"Forward_Scatter",
                                    ifelse(grepl('SSC',name),"Side_Scatter",
                                           ifelse(grepl('node|cluster',name),"FlowSOM_derivative","Unmixed_Fluorescence"))))]
    parms.type<-as.list(stats::setNames(types,nm=paste0("$P",seq(parms.adf[,.N]),'TYPE')))
    volt.names<-paste0("$P",grep("fluorescence|scatter",parms.type,ignore.case = T),'V')
    parms.volts<-as.list(stats::setNames(rep(0,length(volt.names)),nm=volt.names))
    parms.list<-c(parms.list,parms.type,parms.volts)
  }
  return(parms.list)
}
#' @title Convert a `data.table` into a new .fcs file
#' @description
#' After reading in (usually large amounts of) .fcs data and converting to a `data.table`, the process can be reversed and a new .fcs file created; new keywords are written to the header section and parameters updated using `flowCore` functions.
#'
#'
#' @param dt Object returned from `fcs.to.dt`.
#' @param fcs.files.dt Object returned from `get.fcs.file.dt`.
#' @param parse.by Character string; a named column in `dt` that contains unique identifiers; used to both index the `dt` and generate individual, unique .fcs file names.
#' @param type.convert Argument as defined in `parms.adf.dt`
#' @param name.fix Argument as defined in `parms.adf.dt`
#' @param cyto.method Argument as defined in `parms.list.from.adf`
#' @param write Logical; default `TRUE`. A new .fcs file will be written to `out.dir`, using the 'batch' column from `dt` as a sub-directory.
#' @param out.dir File path; default "./data_modified". As files will be modified with both new data columns (`FlowSOM` derived) and keywords/meta-data, they should be written to a "./data_modified" folder so that the source is left intact.
#'
#' @return If `write` is set to `FALSE`, a `flowCore::flowFrame` will be returned; otherwise, a new .fcs file will be written.
#' @export
#'
#'
dt.to.fcs<-function(dt,fcs.files.dt,parse.by='sample.id',type.convert=NULL,name.fix=NULL,cyto.method="spectral",write=T,out.dir="./data_modified"){
  #for R CMD check; data.table vars
  sample.id<-f.path<-batch<-NULL
  #
  drop.terms<-c(
    paste0(c('BEGIN','END'),c('ANALYSIS')),
    paste0(c('BEGIN','END'),c('DATA')),
    paste0(c('BEGINS','ENDS'),c('TEXT')),
    "MODE",
    'BYTEORD',
    'DATATYPE',
    'PAR',
    'TOT',
    'NEXTDATA',
    'DISPLAY'
  )
  #
  parse.ids<-dt[,unique(get(parse.by))]
  #
  for(i in parse.ids){
    parms.adf<-parms.adf.dt(dt[sample.id %in% i],type.convert = type.convert,name.fix = name.fix)
    parms.list<-parms.list.from.adf(parms.adf,cyto.method = cyto.method)
    keywords<-as.list(flowCore::read.FCSheader(fcs.files.dt[sample.id %in% i,f.path])[[1]])
    if(any(grepl("\\$P[0-9]+V",names(keywords)))){
      if(any(unlist(keywords[grep("\\$P[0-9]+V",names(keywords))])!=0)){
        volt.values<-unlist(keywords[grep("\\$P[0-9]+V",names(keywords),value = T)])
        volt.values<-volt.values[volt.values!=0]
        volts.tmp<-stats::setNames(volt.values,nm=sub("_B","B",gsub("-","_",unlist(keywords[sub('V','N',names(volt.values))]))))
        #
        for(v in names(volts.tmp)){
          parms.list[[sub("N","V",names(grep(v,parms.list,value = T)))]]<-volts.tmp[[v]]
        }
      }
    }
    keywords.inherit<-keywords[grep("\\$P[0-9]+",names(keywords),value = T,invert = T)]
    keywords.inherit<-keywords.inherit[grep(paste0(drop.terms,collapse = "|"),names(keywords.inherit),value = T,invert = T)]
    #
    keywords.add<-grep(paste0(sub("\\$","",names(keywords.inherit)),collapse = "|"),names(fcs.files.dt[sample.id %in% i]),value = T,invert = T)
    keywords.add<-grep(paste0(c(drop.terms,'f.path','file.size.MB'),collapse = "|"),keywords.add,invert = T,value = T)
    keywords.inherit<-c(keywords.inherit,as.list(fcs.files.dt[sample.id %in% i,keywords.add,with=F]))
    keywords.inherit[["modified"]]<-'TRUE'
    keywords.inherit[["modified.by"]]<-'SOMnambulate::dt.to.fcs'
    keywords.inherit[["modified.date"]]<-Sys.Date()
    keywords.inherit<-keywords.inherit[which(!is.na(keywords.inherit))]
    #spillover fix
    # spill<-parms.list[names(parms.list) %in% paste0(stringr::str_extract(names(grep("unmixed",parms.list,ignore.case = T,value = T)),"\\$P[0-9]+"),'N')]
    # spill.mat<-matrix(0,nrow=length(spill),ncol=length(spill),dimnames = list(NULL,unlist(spill)))
    # diag(spill.mat)<-1
    # keywords.inherit$`$SPILLOVER`<-spill.mat
    keywords.inherit$`$SPILLOVER`<-NULL
    #
    parms.adf<-data.frame(parms.adf)
    rownames(parms.adf)<-paste0("$P",rownames(parms.adf))
    #
    exprs<-as.matrix(
      dt[sample.id %in% i,.SD,.SDcols = c(dt[,names(.SD),.SDcols = is.numeric],grep(type.convert,names(dt),value = T))][
        ,lapply(.SD,as.numeric)]
    )
    colnames(exprs)<-parms.adf$name
    #
    fcs.from.dt<-methods::new("flowFrame",
                              exprs=exprs,
                              parameters=Biobase::AnnotatedDataFrame(parms.adf),
                              description=parms.list
    )
    flowCore::keyword(fcs.from.dt) <- c(flowCore::keyword(fcs.from.dt),
                                        keywords.inherit)
    flowCore::keyword(fcs.from.dt)[["$FIL"]]<-paste0(flowCore::keyword(fcs.from.dt)[["sample.id"]],".fcs")
    #
    if(!write){
      return(fcs.from.dt)
    }else{
      out.batch<-fcs.files.dt[sample.id %in% i,as.character(batch)]
      out.path<-file.path(out.dir,out.batch)
      out.name<-flowCore::keyword(fcs.from.dt)[["$FIL"]]
      message(
        paste(
          paste("Writing the following .fcs file:", out.name),
          paste("To the following directory:",out.path),
          sep = "\n"
        )
      )
      if(!dir.exists(out.path)){
        dir.create(out.path, recursive = T)
      }
      flowCore::write.FCS(fcs.from.dt, file.path(out.path,out.name))
    }
    #
  }
}

