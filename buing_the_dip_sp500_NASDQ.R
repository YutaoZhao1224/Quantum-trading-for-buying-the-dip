# 加载必要的包
library(quantmod)
library(TTR)
library(dplyr)
library(rvest)


# 1. 获取标普500成分股列表
#url <- "https://en.wikipedia.org/wiki/List_of_S%26P_500_companies"
#sp500_table <- read_html(url) %>% html_nodes("table") %>% .[[1]] %>% html_table()

# 提取股票代码（Ticker）
#sp500_tickers <- sp500_table$Symbol
#sp500_tickers <- gsub("\\.", "-", sp500_tickers)

#sp500 <- data.frame(Symbol = sp500_tickers)
#write.table(sp500, file = "~/Sp500.csv", quote = F, row.names = F)
dir <- "D:/bioinformatics/Quantum Trading/CD_Day Screening/"
dir.create(dir)
BARSLAST <- function(X) {
  last_idx <- rep(NA, length(X))  # 结果向量，初始化为 NA
  last_position <- NA  # 记录上次 TRUE 发生的位置
  
  for (i in seq_along(X)) {
    if (!is.na(X[i]) && X[i]) {
      last_position <- i  # 更新上次 TRUE 发生的位置
    }
    if (!is.na(last_position)) {
      last_idx[i] <- i - last_position  # 计算上次 TRUE 到当前的周期数
    }
  }
  return(last_idx)
}
### 默认五天之内有抄底信号
Choose_stock <- function(stock_name, time_limit = 5, Vol_ratio_period = 5){
  # 3. 设置日期范围
  end_date <- Sys.Date()  # 当前日期
  start_date <- end_date - 200  # 过去200天（包括今天）
  # 获取股票数据
  getSymbols(stock_name, src = "yahoo", from = start_date, to = end_date, auto.assign = TRUE)
  stock_data <- get(stock_name)
  
  # 提取收盘价和量
  close_price <- Cl(stock_data)
  volume <- Vo(stock_data)
  #量比
  len_volume <- length(volume)
  avg_volume <- rep(NA,Vol_ratio_period - 1)
  for (i in Vol_ratio_period:len_volume){
    avg_volume_value <- mean(volume[(i-Vol_ratio_period+1):i])
    avg_volume <- c(avg_volume, avg_volume_value)
  }
  volume_ratio <- volume / (avg_volume+1)
  
  # 计算MACD指标
  #macd <- MACD(close_price, nFast = 12, nSlow = 26, nSig = 9, maType = "EMA")
  #macd <- na.omit(macd)
  DIF <- EMA(close_price[,1],12) - EMA(close_price[,1],26)
  DEA <- EMA(DIF,9)
  MACD_value <- (DIF-DEA) *2
  D <- DIF  # MACD线
  A <- DEA  # 信号线
  M <- (D - A) * 2  # MACD柱状图
  macd <- data.frame(close_price, D,A,M)
  names(macd) <- c("close_price","D","A","M")
  # 计算N1和MM1
  Death_cross <- lag(M,1) >= 0 & M < 0 ### 死叉的时间 TRUE
  Gold_cross <- lag(M,1) <= 0 & M > 0 ### 金叉的时间点 FALSE
  macd <- cbind(macd,BARSLAST(Death_cross), BARSLAST(Gold_cross)) 
  names(macd) <- c("close_price","D","A","M", "N1", "MM1")
  #N1 <- which(diff(sign(M)) == -2)[1]  # 最近一次MACD柱状图由正转负
  #MM1 <- which(diff(sign(M)) == 2)[1]  # 最近一次MACD柱状图由负转正
  
  
  len <- length(DIF)
  # 计算 CC1：最近 N1+1 个周期内的最低价
  last_na_position_N1 <- max(which(is.na(macd$N1)))
  CC1 <- rep(NA,last_na_position_N1)
  for (i in (last_na_position_N1+1):len){
    min_price <- min(macd$close_price[(i-macd$N1[i]-1):i])
    CC1 <- c(CC1,min_price)
  }
  macd$CC1 <- CC1
  # 计算 CC2：CC1 向前延迟 MM1+1 周期
  last_na_position_CC1 <- max(which(is.na(macd$CC1)))
  last_na_position_MM1 <- max(which(is.na(macd$MM1)))
  last_na_position_CC1_MM1 <- max(last_na_position_CC1, last_na_position_MM1)
  CC2 <- rep(NA,last_na_position_CC1_MM1)
  for (i in (last_na_position_CC1_MM1 + 1):len){
    CC2_value <- macd$CC1[i-(macd$MM1[i])-1]
    CC2 <- c(CC2,CC2_value)
  }
  macd$CC2 <- CC2
  # 计算 CC3：CC2 向前延迟 MM1+1 周期
  last_na_position_CC2 <- max(which(is.na(macd$CC2)))
  last_na_position_MM1 <- max(which(is.na(macd$MM1)))
  last_na_position_CC2_MM1 <- max(last_na_position_CC2, last_na_position_MM1)
  CC3 <- rep(NA,last_na_position_CC2_MM1)
  for (i in (last_na_position_CC2_MM1 + 1):len){
    CC3_value <- macd$CC2[i-(macd$MM1[i])-1]
    CC3 <- c(CC3,CC3_value)
  }
  macd$CC3 <- CC3
  
  # 计算价格和MACD极值
  # DIFL1: 最近N1+1个周期内D的最低值
  last_na_position_N1 <- max(which(is.na(macd$N1)))
  last_na_position_D <- max(which(is.na(macd$D)))
  last_na_position_N1_D <- max(last_na_position_N1, last_na_position_D)
  DIFL1 <- rep(NA,last_na_position_N1_D)
  for (i in (last_na_position_N1_D+1):len){
    min_D <- min(macd$D[(i-macd$N1[i]-1):i])
    DIFL1 <- c(DIFL1,min_D)
  }
  macd$DIFL1 <- DIFL1
  # DIFL2:DIFL1在MM1+1周期前的值
  last_na_position_DIFL1 <- max(which(is.na(macd$DIFL1)))
  last_na_position_MM1 <- max(which(is.na(macd$MM1)))
  last_na_position_DIFL1_MM1 <- max(last_na_position_DIFL1, last_na_position_MM1)
  DIFL2 <- rep(NA,last_na_position_DIFL1_MM1)
  for (i in (last_na_position_DIFL1_MM1 + 1):len){
    DIFL2_value <- macd$DIFL1[i-(macd$MM1[i])-1]
    DIFL2 <- c(DIFL2,DIFL2_value)
  }
  macd$DIFL2 <- DIFL2
  # DIFL3:DIFL2在MM1+1周期前的值
  last_na_position_DIFL2 <- max(which(is.na(macd$DIFL2)))
  last_na_position_MM1 <- max(which(is.na(macd$MM1)))
  last_na_position_DIFL2_MM1 <- max(last_na_position_DIFL2, last_na_position_MM1)
  DIFL3 <- rep(NA,last_na_position_DIFL2_MM1)
  for (i in (last_na_position_DIFL2_MM1 + 1):len){
    DIFL3_value <- macd$DIFL2[i-(macd$MM1[i])-1]
    DIFL3 <- c(DIFL3,DIFL3_value)
  }
  macd$DIFL3 <- DIFL3
  ### lagM: M的前一天值
  last_na_position_M <- max(which(is.na(macd$M)))
  last_na_position_M <- max(2,last_na_position_M)
  lagM <- rep(NA,last_na_position_M)
  for (i in (last_na_position_M + 1):len){
    lagM_value <- macd$M[i-1]
    lagM <- c(lagM,lagM_value)
  }
  macd$lagM <- lagM
  
  
  macd <- mutate(macd, AAA = (CC1 < CC2) & (DIFL1 > DIFL2) & (lagM < 0) & (D < 0),
                 BBB = (CC1 < CC3) & (DIFL1 < DIFL2) & (DIFL1 > DIFL3) & lagM < 0 & (D < 0),
                 CCC = (AAA | BBB) & (D < 0))
  ### lagCCC, CCC前一天值 or lag(macd$CCC, 1)
  last_na_position_CCC <- max(which(is.na(macd$CCC)))
  last_na_position_CCC <- max(2,last_na_position_CCC)
  lagCCC <- rep(NA,last_na_position_CCC)
  for (i in (last_na_position_CCC + 1):len){
    lagCCC_value <- macd$CCC[i-1]
    lagCCC <- c(lagCCC,lagCCC_value)
  }
  macd$lagCCC <- lagCCC
  macd <- mutate(macd, LLL = lagCCC & CCC)
  macd$lagAAA <- lag(macd$AAA, 1)
  macd$lagBBB <- lag(macd$BBB, 1)
  macd$lagD <- lag(macd$D, 1)
  macd <- mutate(macd, XXX = (lagAAA & (DIFL1 <= DIFL2) & (D < A)) | (lagBBB & (DIFL1 <= DIFL3) & (D < A)))
  macd <- mutate(macd, JJJ = lagCCC & (abs(lagD) >= (abs(D) * 1.01)))
  macd$lagJJJ <- lag(macd$JJJ, 1)
  macd <- mutate(macd, DXDX = (lagJJJ == F) & JJJ)
  macd$Vol <- volume
  macd$vol_ratio <- volume_ratio
  macd_filtered <- macd[which(macd$DXDX==TRUE),]
  if(nrow(macd_filtered)>0){
    DayTime <- tail(macd_filtered,1)%>% row.names()
    DayDate <- strsplit(DayTime, " ")[[1]][1]
    date_to_check <- as.Date(DayDate)
    # 获取当前日期
    current_date <- Sys.Date()
    
    # 判断是否在5天内
    is_within_5_days <- (current_date - date_to_check) <= time_limit
    if(is_within_5_days){
      print(paste0(stock_name, ":chosen"))
      VOL <- tail(macd_filtered,1) %>% pull(Vol)
      VOL <- VOL/1000000
      VOL_ratio <- tail(macd_filtered,1) %>% pull(vol_ratio)
      price <- tail(macd_filtered,1) %>% pull(close_price)
      return(c(DayDate, stock_name, price, VOL, VOL_ratio))
    }  else {
      print(paste0(stock_name, ":not chosen"))
      return(NULL)
    }
  } else {
    print(paste0(stock_name, ":not chosen"))
    return(NULL)
  }
}

Sp500 <- read.csv("D:/bioinformatics/Quantum Trading/Sp500.csv", sep="") 
sp500_tickers <- Sp500$Symbol
stock_list_sp500 <- c()
for (i in sp500_tickers){
  try({
    print(which(sp500_tickers==i))
    stock_chosen <- Choose_stock(i,5)
    stock_list_sp500<- c(stock_list_sp500,stock_chosen)}, silent = TRUE)
}
sp500_stock_data <- data.frame(
  Date = as.Date(stock_list_sp500[seq(1, length(stock_list_sp500), by = 5)]), # 抄底时间
  Symbol = stock_list_sp500[seq(2, length(stock_list_sp500), by = 5)],  # 股票代码
  Price = as.numeric(stock_list_sp500[seq(3, length(stock_list_sp500), by = 5)]),  # 价格
  Volume = as.numeric(stock_list_sp500[seq(4, length(stock_list_sp500), by = 5)]),  # 交易量
  Volume_Ratio = as.numeric(stock_list_sp500[seq(5, length(stock_list_sp500), by = 5)])  # 量比
)
sp500_stock_data <- sp500_stock_data %>% arrange(desc(Date)) 
write.table(sp500_stock_data, file = paste0(dir, "sp500_stock_data--", sp500_stock_data$Date[1], "--.txt"),
            quote = FALSE, row.names = FALSE, col.names = T, sep = "\t")

NASDQ <- read.csv("D:/bioinformatics/Quantum Trading/NASDQ.csv")
NASDQ_tickers <- NASDQ$Symbol
stock_list_NASDQ <- c()
for (i in NASDQ_tickers){
  try({
    print(which(NASDQ_tickers==i))
    stock_chosen <- Choose_stock(i,5)
    stock_list_NASDQ <- c(stock_list_NASDQ,stock_chosen)}, silent = TRUE)
}
NASDQ_stock_data <- data.frame(
  Date = as.Date(stock_list_NASDQ[seq(1, length(stock_list_NASDQ), by = 5)]), # 抄底时间
  Symbol = stock_list_NASDQ[seq(2, length(stock_list_NASDQ), by = 5)],  # 股票代码
  Price = as.numeric(stock_list_NASDQ[seq(3, length(stock_list_NASDQ), by = 5)]),  # 价格
  Volume = as.numeric(stock_list_NASDQ[seq(4, length(stock_list_NASDQ), by = 5)]),  # 交易量
  Volume_Ratio = as.numeric(stock_list_NASDQ[seq(5, length(stock_list_NASDQ), by = 5)])  # 量比
)
NASDQ_stock_data <- NASDQ_stock_data %>% arrange(desc(Date)) 
write.table(NASDQ_stock_data, file = paste0(dir, "NASDQ_stock_data--", NASDQ_stock_data$Date[1], "--.txt"),
            quote = FALSE, row.names = FALSE, col.names = T, sep = "\t")
