#!/usr/bin/env Rscript

# 1. ЗАГРУЗКА БИБЛИОТЕК
library(jsonlite)
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(VariantAnnotation)
library(httr)
library(tidyr)
library(purrr)

# BiocManager::install(c("EnsDb.Hsapiens.v86", "ensembldb"))
library(EnsDb.Hsapiens.v86)
library(ensembldb)


# Создаем объект, который искал код
hg38_genome <- BSgenome.Hsapiens.UCSC.hg38
# Принудительно устанавливаем стиль UCSC (chr1, chr2...), так как BSgenome его использует
seqlevelsStyle(hg38_genome) <- "UCSC" 

# --- ФУНКЦИЯ РАСЧЕТА HGVS.p (Улучшенная версия для интеграции) ---
calculate_hgvs_p_integrated <- function(chrom, pos, ref, alt) {
  aa_map <- c("M"="Met","T"="Thr","A"="Ala","R"="Arg","N"="Asn","D"="Asp","C"="Cys",
              "Q"="Gln","E"="Glu","G"="Gly","H"="His","I"="Ile","L"="Leu","K"="Lys",
              "F"="Phe","P"="Pro","S"="Ser","W"="Trp","Y"="Tyr","V"="Val")
  
  # Форматирование хромосомы для Ensembl (без 'chr')
  chrom_ens <- gsub("chr", "", as.character(chrom))
  chrom_ucsc <- if(!grepl("chr", chrom)) paste0("chr", chrom) else chrom
  
  # Поиск позиции в белке
  gn_range_ens <- GRanges(chrom_ens, IRanges(start = pos, width = nchar(ref)))
  # Подавляем предупреждения о несоответствии стилей хромосом
  pr_pos <- suppressWarnings(genomeToProtein(gn_range_ens, EnsDb.Hsapiens.v86))
  res_df <- as.data.frame(unlist(pr_pos))
  
  if (nrow(res_df) == 0) return(NA_character_)
  
  # Берем первый транскрипт (или можно сопоставить с tx_id из vcf)
  res_df <- res_df[1, ] 
  
  # Работа с кодонами через BSgenome
  # Берем окрестность, чтобы захватить кодон
  codon_range <- GRanges(chrom_ucsc, IRanges(start = pos - 1, end = pos + 1))
  ref_codon <- getSeq(BSgenome.Hsapiens.UCSC.hg38, codon_range)
  
  alt_codon_seq <- as.character(ref_codon)
  # Заменяем центральный нуклеотид (упрощенная логика для SNP)
  if(nchar(alt) == 1 && nchar(ref) == 1) {
    substr(alt_codon_seq, 2, 2) <- alt
    alt_codon <- DNAString(alt_codon_seq)
    
    ref_aa <- as.character(translate(ref_codon))
    alt_aa <- as.character(translate(alt_codon))
    
    # Сборка HGVS.p
    ref_3 <- ifelse(!is.na(aa_map[ref_aa]), aa_map[ref_aa], ref_aa)
    alt_3 <- ifelse(!is.na(aa_map[alt_aa]), aa_map[alt_aa], alt_aa)
    return(paste0("p.", ref_3, res_df$start, alt_3))
  } else {
    return("p.ComplexIndel") # Для инваров логика сложнее
  }
}

# *************************************************************************************
# *************************** 1. ЗАГРУЗКА КОНФИГУРАЦИИ ********************************
# *************************************************************************************

message("Reading configuration settings...")

# Загрузка конфига (путь можно изменить на "./config.json", если файл в текущей папке)
config_path <- "~/R/bind/GDF2425Y/scripts/config.json"
config <- fromJSON(config_path)

# Проверка наличия нужных полей в конфигурации
if (!all(c("project_dir", "nuclear_config") %in% names(config))) {
  stop("Error: Configuration file missing required fields!")
}

config$nuclear_config$dir_path <- sub("${project_dir}", config$project_dir, config$nuclear_config$dir_path, fixed = TRUE)
config$nuclear_config$output_dir <- sub("${project_dir}", config$project_dir, config$nuclear_config$output_dir, fixed = TRUE)

# Другие настройки
genes_symbol <- config$nuclear_config$genes_symbol
max_af_nfe <- config$nuclear_config$max_af_nfe
consequences <- config$nuclear_config$consequences

# Создаем директорию для результатов
dir.create(config$nuclear_config$output_dir, showWarnings=FALSE, recursive=TRUE)

# Выводим текущие настройки
message("Configuration:")
message("Project Directory:", config$project_dir)
message("Samples Directory:", config$samples_dir)
message("Reference Genome File:", config$refgen_file)
message("Rename Chromosome Map:", config$rename_chromosomes_map)
message("Number of Threads BWA:", config$num_threads_bwa)
message("Number of Threads Sort:", config$num_threads_sort)
message("Number of Threads GATK:", config$num_threads_gatk)
message("Log to Terminal:", config$log_to_terminal)

# Ядерные настройки
message("Nuclear Config:")
message("Genes Symbol:", paste(config$nuclear_config$genes_symbol, collapse = ", "))
message("Max AF NFE:", config$nuclear_config$max_af_nfe)
message("Consequences:", paste(config$nuclear_config$consequences, collapse = ", "))
message("Directory Path:", config$nuclear_config$dir_path)
message("Output Directory:", config$nuclear_config$output_dir)

###########################################################################################
# Начало основного процесса
message("[", Sys.time(), "] Starting annotation pipeline")

# *************************************************************************************
# *************************** 2. ПОДГОТОВКА СПРАВОЧНОЙ ИНФОРМАЦИИ (HG38) **************
# *************************************************************************************


# message("Preparing CLINVAR...")
# 
# # Указываем точный путь к скачанным файлам
# clinvar_path <- "~/R/bind/GDF2425Y/refgen/clivar/clinvar_20260120.vcf.gz"

message("Preparing CLINVAR...")
clinvar_path <- "/home/rstudio/R/bind/GDF2425Y/refgen/clivar/clinvar_20260120.vcf.gz"

if (!file.exists(paste0(clinvar_path, ".tbi"))) {
  # Если индекса нет, лучше остановить скрипт, иначе VariantAnnotation упадет позже
  stop("Ошибка: Файл индекса .tbi не найден по пути: ", clinvar_path)
} 

message("Preparing gene models and transcript ranges...")
# 
# message("Preparing CLINVAR...")
# clinvar_path <- "~/R/bind/GDF2425Y/refgen/clivar/clinvar_20260120.vcf.gz"
# 
# if (!file.exists(paste0(clinvar_path, ".tbi"))) {
#   stop("КРИТИЧЕСКАЯ ОШИБКА: Файл индекса .tbi не найден! Без него поиск невозможен.")
# } # <--- 
# Проверка связи с индексом (tbi должен лежать рядом с тем же именем)

# if (!file.exists(paste0(clinvar_path, ".tbi"))) {
#   message("Внимание: Файл индекса .tbi не найден в /refgen/clivar/. Проверьте имя!")

# message("Preparing gene models and transcript ranges...")

 
  
# Получение Entez ID для символов генов
genes_id <- AnnotationDbi::select(org.Hs.eg.db, 
                   keys = genes_symbol, 
                   keytype = "SYMBOL", 
                   columns = c("ENTREZID", "SYMBOL"))

# Использование референсной базы HG38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Создание списка транскриптов, сгруппированных по генам
genes_txGRL <- transcriptsBy(txdb, "gene")[genes_id$ENTREZID]

# Генерация плоской таблицы транскриптов для быстрого сопоставления
genes_tx_list <- lapply(names(genes_txGRL), function(id) {
  df <- as.data.frame(genes_txGRL[[id]])
  df$GENE_SYMBOL <- genes_id$SYMBOL[match(id, genes_id$ENTREZID)]
  return(df)
})
genes_txDF <- do.call(rbind, genes_tx_list)

# Создание GRanges объекта только для целевых транскриптов
genes_txGR <- transcripts(txdb, columns = c("tx_id", "tx_name"), 
                          filter = list(tx_name = genes_txDF$tx_name), 
                          use.names = TRUE)

# *************************************************************************************
# *************************** 3. СБОР ДАННЫХ GNOMAD V4  *******************************
# *************************************************************************************

message("Fetching data from gnomAD database...")

# Функция для получения данных о вариантах генов из gnomAD
fetch_gene_variants <- function(gene_symbol) {
  query <- paste0('query VariantsInGene {',
                  '  gene(gene_symbol: "', gene_symbol, '", reference_genome: GRCh38) {',
                  '    variants(dataset: gnomad_r4) {',
                  '      variant_id, chrom, pos, ref, alt, rsids, transcript_id, transcript_version, hgvs, hgvsc, hgvsp, consequence, flags, ',
                  '      joint { ac, an, populations { id, ac, an, homozygote_count, hemizygote_count } }',
                  '    }',
                  '  }',
                  '}')
  
  url <- "https://gnomad.broadinstitute.org/api/"
  headers <- add_headers(Content_Type = "application/json")
  body <- list(query = query)
  
  response <- POST(url, body = body, encode = "json", .headers = headers)
  
  if (response$status_code != 200) {
    stop("Ошибка при отправке запроса:", content(response, "text"))
  }
  
  raw_data <- content(response, "text")
  result <- fromJSON(raw_data)
  df <- as.data.frame(result$data$gene$variants, stringsAsFactors = FALSE)
  df$gene_symbol <- gene_symbol
  return(df)
}

# Собираем данные по всем генам
gene_var_gad_list <- lapply(genes_symbol, fetch_gene_variants)
gene_var_gadDF    <- bind_rows(gene_var_gad_list)

# Расчёт общей частоты аллелей (Joint AF)
gene_var_gadDF$joint_af <- with(gene_var_gadDF, ifelse(!is.na(joint$ac) & !is.na(joint$an) & joint$an > 0, joint$ac / joint$an, NA))

# Извлечение статистики по популяции NFE (European Non-Finnish)
message("Processing NFE population statistics (Fast mode)...")
# 1. Извлекаем все данные популяций и присваиваем им имена из variant_id
pops_list <- gene_var_gadDF$joint$populations
names(pops_list) <- gene_var_gadDF$variant_id

# 2. Быстро объединяем в одну таблицу (это намного быстрее lapply)
nfe_pop_data <- dplyr::bind_rows(pops_list, .id = "variant_id") %>%
  # Фильтруем только нужную популяцию
  dplyr::filter(id == "nfe") %>%
  # Переименовываем, чтобы не было конфликтов с основными колонками ac/an
  dplyr::select(variant_id, 
                ac_nfe = ac, 
                an_nfe = an, 
                hom_nfe = homozygote_count)

# 3. Соединяем с основной таблицей gnomAD
gnomad_ref <- gene_var_gadDF %>%
  dplyr::left_join(nfe_pop_data, by = "variant_id") %>%
  dplyr::mutate(
    # Общая частота (Joint AF) из верхнего уровня
    joint_af = if_else(!is.na(joint$an) & joint$an > 0, joint$ac / joint$an, as.numeric(NA)),
    # Частота только для NFE
    af_nfe = if_else(!is.na(an_nfe) & an_nfe > 0, ac_nfe / an_nfe, 0)
  ) %>%
  # Оставляем финальный набор колонок
  dplyr::select(
    variant_id, rsids, transcript_id, transcript_version, hgvs, consequence, 
    gene_symbol, joint_af, af_nfe, ac_nfe, an_nfe, hom_nfe
  ) %>%
  # Убираем дубликаты, если они возникли
  dplyr::distinct(variant_id, .keep_all = TRUE)

# *************************************************************************************
# *************************** 4. ЦИКЛ ОБРАБОТКИ VCF-ФАЙЛОВ ****************************
# *************************************************************************************

message("Starting processing of VCF files...")

vcf_files <- list.files(path = config$nuclear_config$dir_path, pattern = "\\.vcf$", full.names = TRUE)

for(file in vcf_files) {
  sample_base <- tools::file_path_sans_ext(basename(file))
  message("\n[", Sys.time(), "] Processing: ", sample_base)
  
  # 1. Сначала определяем и создаем директорию
  sample_dir <- file.path(config$nuclear_config$output_dir, sample_base)
  dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 2. Теперь создаем базовое имя для файлов внутри этой папки
  file_out_base <- file.path(sample_dir, sample_base) 
  
  # ... далее чтение VCF и остальной код ...
  # sample_base <- tools::file_path_sans_ext(basename(file))
  # message("\n[", Sys.time(), "] Processing: ", sample_base)
  # 
  # # --- СОЗДАНИЕ ПЕРСОНАЛЬНОЙ ПАПКИ ОБРАЗЦА ---
  # sample_dir <- file.path(config$nuclear_config$output_dir, sample_base)
  # # ДОБАВЬТЕ ЭТУ СТРОКУ СЮДА:
  # file_out_base <- file.path(sample_dir, sample_base) 
  # 
  # message("\n[", Sys.time(), "] Processing: ", sample_base)
  # # ... остальной код цикла ...
  # dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Чтение VCF-файла
  vcf_raw <- VariantAnnotation::readVcf(file, genome = "hg38")
  
  # Фильтрация по пересечениям с транскриптами
  vcf_by_tx_collapsed <- subsetByOverlaps(vcf_raw, genes_txGR)
  if (nrow(vcf_by_tx_collapsed) == 0) {
    message("Warning: No variants found in target genes for ", sample_base, ". Skipping.")
    next
  }
  
  # Развертывание мультивалентных вариантов
  vcf_by_tx <- VariantAnnotation::expand(vcf_by_tx_collapsed) 
  vcf_df <- as.data.frame(SummarizedExperiment::rowRanges(vcf_by_tx))
  
  # --- НОВЫЙ СПОСОБ ИЗВЛЕЧЕНИЯ ГЕНОТИПОВ (LONG FORMAT) ---
  matGT <- VariantAnnotation::geno(vcf_by_tx)$GT
  matDP <- if ("DP" %in% names(VariantAnnotation::geno(vcf_by_tx))) VariantAnnotation::geno(vcf_by_tx)$DP else NULL
  
  # Вместо цикла по колонкам, создаем одну запись для каждого образца
  # Если образец в VCF один, это просто создаст колонки GT и DP
  if (ncol(matGT) >= 1) {
    # Берем первый образец (в вашем случае это и есть целевой образец из файла)
    vcf_df$GT <- matGT[, 1]
    vcf_df$DP <- if (!is.null(matDP)) matDP[, 1] else NA
    
    vcf_df$GT_Desc <- dplyr::case_when(
      vcf_df$GT %in% c("0/1", "1/0", "0|1", "1|0") ~ "Heterozygous",
      vcf_df$GT %in% c("1/1", "1|1") ~ "Homozygous Mutant",
      vcf_df$GT %in% c("0/0", "0|0") ~ "Wild Type",
      TRUE ~ "Other"
    )
  }
  
  # Удаляем старую логику, которая создавала колонки типа ASA250720_GT
  
  
  # Аннотация транскриптов и генов
  hits <- findOverlaps(SummarizedExperiment::rowRanges(vcf_by_tx), genes_txGR)
  vcf_df$tx_name <- NA_character_
  vcf_df$GENE_SYMBOL <- NA_character_
  
  if(length(hits) > 0) {
    vcf_df$tx_name[queryHits(hits)] <- genes_txGR$tx_name[subjectHits(hits)]
    vcf_df$GENE_SYMBOL[queryHits(hits)] <- genes_txDF$GENE_SYMBOL[match(vcf_df$tx_name[queryHits(hits)], genes_txDF$tx_name)]
  }
  
  # Создание варианта ID для стыковки с gnomAD
  vcf_df$variant_id <- paste(gsub("chr", "", as.character(vcf_df$seqnames)), 
                             vcf_df$start, vcf_df$REF, as.character(vcf_df$ALT), sep = "-")
  
  # Объединение с данными gnomAD (используем dplyr явно)
  vcf_annotated <- vcf_df %>% dplyr::left_join(gnomad_ref, by = "variant_id")
  
  # Дополнение информацией
  vcf_annotated <- vcf_annotated %>%
    dplyr::mutate(
      transcript_full = paste(transcript_id, transcript_version, sep = "."),
      af_nfe = dplyr::if_else(!is.na(an_nfe) & an_nfe > 0, ac_nfe / an_nfe, 0)
    )
  
  # Очистка списковых колонок
  list_cols <- names(vcf_annotated)[sapply(vcf_annotated, is.list)]
  for(col in list_cols) {
    vcf_annotated[[col]] <- sapply(vcf_annotated[[col]], function(x) {
      if(is.null(x) || length(x) == 0) return(NA_character_)
      paste(as.character(unlist(x)), collapse=",")
    })
  }
  
  # Фильтрация приоритетных вариантов
  priority_df <- vcf_annotated %>%
    dplyr::filter(af_nfe < config$nuclear_config$max_af_nfe | is.na(af_nfe)) %>%
    dplyr::filter(consequence %in% config$nuclear_config$consequences)
  

  if (nrow(priority_df) > 0) {
    message("Predicting pathogenicity for priority variants...")
    
    # 1. Подготовка GRanges и синхронизация (стандарт UCSC для hg38)
    priority_gr <- GenomicRanges::makeGRangesFromDataFrame(priority_df, keep.extra.columns = TRUE)
    priority_gr <- GenomeInfoDb::keepStandardChromosomes(priority_gr, pruning.mode = "coarse")
    GenomeInfoDb::seqlevelsStyle(priority_gr) <- "UCSC"
    
    txdb_sync <- GenomeInfoDb::keepStandardChromosomes(txdb, pruning.mode = "coarse")
    GenomeInfoDb::seqlevelsStyle(txdb_sync) <- "UCSC"
    GenomeInfoDb::seqlevels(priority_gr) <- GenomeInfoDb::seqlevels(txdb_sync)
    
    # 2. Расчет аминокислотных замен
    coding_preds <- tryCatch({
      VariantAnnotation::predictCoding(
        query = priority_gr, 
        subject = txdb_sync, 
        seqSource = hg38_genome, 
        varAllele = Biostrings::DNAStringSet(priority_df$ALT)
      )
    }, error = function(e) {
      message("Warning: predictCoding failed: ", e$message)
      return(NULL)
    })
    
    # 3. Интеграция результатов
    if (!is.null(coding_preds) && length(coding_preds) > 0) {
      coding_df <- as.data.frame(coding_preds) %>%
        dplyr::mutate(
          seqnames = gsub("chr", "", as.character(seqnames)),
          # Извлекаем позицию аминокислоты (берем только первое значение)
          # p_pos = sapply(PROTEINLOC, function(x) if(length(x) > 0) head(BiocGenerics::start(x), 1) else NA)
          p_pos = sapply(PROTEINLOC, function(x) if(length(x) > 0) as.integer(unlist(x)[1]) else NA)
          
        ) %>%
        dplyr::select(
          seqnames, start, 
          AA_CALC = REFAA, 
          VAR_CALC = VARAA,
          POS_CALC = p_pos,
          CONS_TYPE_CALC = CONSEQUENCE
        ) %>%
        dplyr::distinct(seqnames, start, .keep_all = TRUE)
      
      # ПРОВЕРКА И СОЗДАНИЕ КОЛОНОК (чтобы coalesce не падал)
      for(col in c("AA_CHANGE", "VAR_AA", "PROTEIN_POS", "CONSEQUENCE_TYPE")) {
        if(!(col %in% colnames(priority_df))) priority_df[[col]] <- NA_character_
      }
      
      priority_df <- priority_df %>%
        dplyr::mutate(seqnames = gsub("chr", "", as.character(seqnames))) %>%
        dplyr::left_join(coding_df, by = c("seqnames", "start")) %>%
        dplyr::mutate(
          AA_CHANGE = dplyr::coalesce(as.character(AA_CHANGE), as.character(AA_CALC)),
          VAR_AA = dplyr::coalesce(as.character(VAR_AA), as.character(VAR_CALC)),
          PROTEIN_POS = dplyr::coalesce(as.character(PROTEIN_POS), as.character(POS_CALC)),
          CONSEQUENCE_TYPE = dplyr::coalesce(as.character(CONSEQUENCE_TYPE), as.character(CONS_TYPE_CALC))
        ) %>%
        dplyr::select(-dplyr::ends_with("_CALC")) # Удаляем временные колонки
    }
    
    # 4. ФИНАЛЬНАЯ СБОРКА И ПРИОРИТИЗАЦИЯ
    priority_df <- priority_df %>%
      dplyr::mutate(
        HGVS_P = dplyr::case_when(
          !is.na(hgvs) & hgvs != "" ~ as.character(hgvs),
          !is.na(AA_CHANGE) & AA_CHANGE != "NA" ~ paste0("p.", AA_CHANGE, PROTEIN_POS, VAR_AA),
          TRUE ~ NA_character_
        ),
        HGVS_P = dplyr::if_else(grepl("NA", HGVS_P), NA_character_, HGVS_P)
      ) %>%
      dplyr::arrange(dplyr::desc(is.na(af_nfe)), af_nfe) %>%
      dplyr::relocate(dplyr::any_of(c("GENE_SYMBOL", "HGVS_P", "GT_Desc", "af_nfe", "consequence")))
  }

  ### CLINVAR BEG
  ### CLINVAR BEG
  if (nrow(priority_df) > 0) {
    message("Consulting ClinVar for priority variants (Advanced Info)...")
    
    priority_gr_cv <- GenomicRanges::makeGRangesFromDataFrame(priority_df)
    GenomeInfoDb::seqlevelsStyle(priority_gr_cv) <- "NCBI" 
    
    # Запрос расширенных полей
    cv_param <- VariantAnnotation::ScanVcfParam(
      which = priority_gr_cv, 
      info = c("CLNSIG", "CLNDN", "CLNREVSTAT", "CLNSIGCONF", "CLNDISDB")
    )
    
    cv_results <- tryCatch({
      VariantAnnotation::readVcf(clinvar_path, "hg38", cv_param)
    }, error = function(e) {
      message("ClinVar lookup skipped: ", e$message)
      return(NULL)
    })
    
    if (!is.null(cv_results) && nrow(cv_results) > 0) {
      cv_info <- VariantAnnotation::info(cv_results)
      cv_ranges <- SummarizedExperiment::rowRanges(cv_results)
      
      cv_df <- data.frame(
        seqnames = gsub("chr", "", as.character(seqnames(cv_ranges))),
        start = start(cv_ranges),
        REF = as.character(ref(cv_results)),
        ALT = as.character(unlist(alt(cv_results))),
        CLINVAR_SIG = sapply(cv_info$CLNSIG, function(x) paste(x, collapse = "/")),
        CLINVAR_DISEASE = sapply(cv_info$CLNDN, function(x) paste(x, collapse = "/")),
        CLINVAR_REV = sapply(cv_info$CLNREVSTAT, function(x) paste(x, collapse = "/")),
        CLINVAR_CONF = sapply(cv_info$CLNSIGCONF, function(x) paste(x, collapse = "/")),
        CLINVAR_DB = sapply(cv_info$CLNDISDB, function(x) paste(x, collapse = "/")),
        stringsAsFactors = FALSE
      ) %>% dplyr::distinct(seqnames, start, REF, ALT, .keep_all = TRUE)
      
      priority_df <- priority_df %>%
        dplyr::left_join(cv_df, by = c("seqnames", "start", "REF", "ALT"))
    } else {
      # ПРАВКА: Создаем ВСЕ колонки, чтобы Rmd не ругался
      priority_df$CLINVAR_SIG <- "not_in_clinvar"
      priority_df$CLINVAR_DISEASE <- NA_character_
      priority_df$CLINVAR_REV <- NA_character_
      priority_df$CLINVAR_CONF <- NA_character_
      priority_df$CLINVAR_DB <- NA_character_
    }
    
    # Выносим всё важное в начало
    priority_df <- priority_df %>%
      dplyr::relocate(dplyr::any_of(c("GENE_SYMBOL", "HGVS_P", "CLINVAR_SIG", 
                                      "CLINVAR_DISEASE", "CLINVAR_REV", "af_nfe")))
  }
  ### CLINVAR END  
    
  # # --- ОЧИСТКА ПРИОРИТЕТНОЙ ТАБЛИЦЫ ДЛЯ МЕДИЦИНСКОГО ОТЧЕТА ---
  # --- ФОРМИРОВАНИЕ ТАБЛИЦЫ ДЛЯ ПЕЧАТИ (ОТЧЕТА) ---
  if (nrow(priority_df) > 0) {
    # Создаем облегченную версию для CSV
    priority_report_df <- priority_df %>%
      dplyr::select(
        # 1. Главное (Клиническая интерпретация)
        GENE_SYMBOL, HGVS_P, CLINVAR_SIG, af_nfe, consequence, GT_Desc,
        
        # 2. Технические данные верификации
        seqnames, start, REF, ALT, QUAL, DP, 
        
        # 3. Идентификаторы
        rsids, transcript_full
      )
    
    # Сохраняем "чистый" CSV для человека
    write.csv(priority_report_df, paste0(file_out_base, "_priority_REPORT.csv"), row.names=FALSE, na="")
  }
  
  # В RDS бандл (для системы) сохраняем ПОЛНЫЙ priority_df со всеми 30+ колонками
  analysis_bundle <- list(
    full_annotated = vcf_annotated,
    priority_variants = priority_df, # Тут все данные сохранены
    sample_id = sample_base,
    analysis_date = Sys.time()
  )
  


# --- СОХРАНЕНИЕ РЕЗУЛЬТАТОВ В ПАПКУ ОБРАЗЦА ---
# Определяем базовое имя внутри папки
  file_out_base <- file.path(sample_dir, sample_base)
  
  # Функция для сохранения во всех форматах
  save_results <- function(df, suffix) {
    path_no_ext <- paste0(file_out_base, suffix)
    # TSV
    write.table(df, paste0(path_no_ext, ".tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")
    # CSV
    write.csv(df, paste0(path_no_ext, ".csv"), row.names=FALSE, na="")
    # RDS
    saveRDS(df, paste0(path_no_ext, ".rds"))
  }
  
  save_results(vcf_annotated, "_full")
  save_results(priority_df, "_priority")

########
  # --- ОБНОВЛЕННОЕ СОХРАНЕНИЕ РЕЗУЛЬТАТОВ (BUNDLE) ---
  file_out_base <- file.path(sample_dir, sample_base)
  
  # 1. Формируем единый пакет данных (Bundle)
  # Здесь мы сохраняем ВСЁ: и полные таблицы, и отфильтрованные
  analysis_bundle <- list(
    full_annotated = vcf_annotated, # Огромная таблица со всеми данными
    priority_variants = priority_df, # Отфильтрованные данные
    sample_id = sample_base,
    analysis_date = Sys.time(),
    settings = config$nuclear_config # Или mito_config
  )
  
  # 2. Сохраняем RDS (это будет основной источник для вашего отчета)
  saveRDS(analysis_bundle, paste0(file_out_base, "_results_bundle.rds"))
  
  # 3. (Опционально) Сохраняем CSV/TSV для быстрой проверки человеком
  # Но теперь отчет Rmd будет брать данные ТОЛЬКО из RDS выше
  write.csv(priority_df, paste0(file_out_base, "_priority.csv"), row.names=FALSE, na="")
  
  
########
  
  message("Processed successfully. Files saved in: ", sample_dir)
  message("Priority variants count: ", nrow(priority_df))
}

message("\n[", Sys.time(), "] Nuclear Pipeline Completed.")

# Скрипт сборки сводной таблицы (Summary Table)
# Запустите этот код после завершения основного цикла обработки.
# 
# library(dplyr)
# library(purrr)

message("\n[", Sys.time(), "] Сборка сводной таблицы по всем образцам...")

# 1. Поиск всех файлов _priority.rds, которые создал основной цикл
summary_files <- list.files(
  path = config$nuclear_config$output_dir, 
  pattern = "_priority\\.rds$", 
  recursive = TRUE, 
  full.names = TRUE
)

# 2. Чтение и объединение
all_priorities <- summary_files %>%
  map_df(function(file_path) {
    # Читаем данные
    df <- readRDS(file_path)
    
    # Извлекаем ID образца из имени файла или папки
    s_id <- basename(dirname(file_path))
    
    # Добавляем колонку с ID образца в самое начало
    if (nrow(df) > 0) {
      df <- df %>% mutate(Sample_ID = s_id) %>% relocate(Sample_ID)
    }
    return(df)
  })

# Оставляем только заполненные и важные колонки в Summary
# all_priorities <- all_priorities %>%
#   dplyr::select(
#     Sample_ID, GENE_SYMBOL, HGVS_P, GT_Desc, af_nfe, consequence,
#     seqnames, start, REF, ALT, QUAL, DP, rsids, transcript_full
#   ) %>%
#   # Сортировка: Сначала по ID пациента, потом по редкости мутации
#   dplyr::arrange(Sample_ID, af_nfe)
# 

# ... (начало блока без изменений: поиск файлов и map_df) ...

# Оставляем только заполненные и важные колонки в Summary
all_priorities <- all_priorities %>%
  dplyr::select(
    Sample_ID, 
    GENE_SYMBOL, 
    HGVS_P, 
    CLINVAR_SIG, 
    CLINVAR_DISEASE,  # Добавлено
    CLINVAR_REV,      # Добавлено (звезды)
    CLINVAR_CONF,     # Добавлено (детали конфликтов)
    CLINVAR_DB,       # Добавлено (OMIM/Orphanet)
    af_nfe, 
    consequence,
    GT_Desc,          # Добавлено обратно для отчета
    seqnames, 
    start, 
    REF, 
    ALT, 
    QUAL, 
    DP, 
    rsids, 
    transcript_full
  ) %>%
  # Сортировка: сначала по гену, затем по редкости мутации
  dplyr::arrange(GENE_SYMBOL, af_nfe)

# 3. Сохранение финальной сводки
if (nrow(all_priorities) > 0) {
  output_summary <- file.path(config$nuclear_config$output_dir, "PROJECT_SUMMARY_PRIORITY_2026.csv")
  
  write.csv(all_priorities, output_summary, row.names = FALSE, na = "")
  
  message("Сводная таблица готова! Поля ClinVar (Advanced) включены.")
  message("Найдено вариантов: ", nrow(all_priorities))
  message("Путь к файлу: ", output_summary)
} else {
  message("Приоритетных вариантов не найдено ни в одном образце.")
}

#----OLD
# all_priorities <- all_priorities %>%
#   dplyr::select(
#     Sample_ID, GENE_SYMBOL, HGVS_P, CLINVAR_SIG, af_nfe, consequence,
#     seqnames, start, REF, ALT, QUAL, DP, rsids, transcript_full
#   ) %>%
#   dplyr::arrange(Sample_ID, af_nfe)
# 
# # 3. Сохранение финальной сводки
# if (nrow(all_priorities) > 0) {
#   output_summary <- file.path(config$nuclear_config$output_dir, "PROJECT_SUMMARY_PRIORITY_2026.csv")
#   
#   # Сортировка: сначала по гену, затем по редкости
#   all_priorities <- all_priorities %>%
#     arrange(GENE_SYMBOL, af_nfe)
#   
#   write.csv(all_priorities, output_summary, row.names = FALSE, na = "")
#   
#   message("Сводная таблица готова! Найдено вариантов: ", nrow(all_priorities))
#   message("Путь к файлу: ", output_summary)
# } else {
#   message("Приоритетных вариантов не найдено ни в одном образце.")
# }

#---END OLD



