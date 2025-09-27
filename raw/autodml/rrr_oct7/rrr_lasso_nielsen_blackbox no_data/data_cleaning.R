# read raw file

df<-read.csv2("monthly_pr_index_lag_imputed_data.csv", sep=",")
#df_old<-read.csv2("final_lag_and_imputed_prices_data.csv", sep=",")

df <- as.data.frame(sapply(df, as.numeric))
summary(df)

#df_old <- as.data.frame(sapply(df_old, as.numeric))
#summary(df_old)

# keep relevant cols only
names(df)

# QUESTION: IS TOTAL SPENT IN SAME SCALE AS BEFORE?
# PRICE INDEXES ARE LOG OF WHAT THEY WERE BEFORE
# MILK_EXP IS WHAT IT WAS BEFORE

df<-dplyr::select(df,
                  #panel_year, #debug
                  month_index, #debug
                  household_code, #i
                  total_spent, #needed to construct LHS vars
                  milk_exp, #as before
                  soda_exp,
                  bread_pr_index, #needed to construct RHS vars
                  butter_pr_index, #log compared to previous dataset; deleted log transform below
                  cereal_pr_index,
                  chips_pr_index,
                  coffee_pr_index,
                  cookies_pr_index,
                  eggs_pr_index,
                  ic_pr_index,
                  milk_pr_index,
                  oj_pr_index,
                  salad_pr_index,
                  soda_pr_index,           
                  soup_pr_index,
                  water_pr_index,
                  yogurt_pr_index) #debug

df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))))

summary(df)

#dropping dropping observations with milk_exp=NA and soda_exp=NA
#df<-df[complete.cases(df),] # so many NAs!

#recoding observations with milk_exp=NA and soda_exp=NA as milk_exp=0 and soda_exp=0
df[is.na(df)] <- 0

# t=1...T
df<-add_count(df,household_code)
df<-df %>% group_by(household_code) %>% mutate(t = row_number())
df<-ungroup(df)

df<-transmute(df,
              hh=as.numeric(as.factor(household_code)), #i
              t=t,
              n=n,
              
              milk_s=milk_exp/total_spent, #LHS vars (shares), like Ying's code
              soda_s=soda_exp/total_spent,
              
              expend=log(total_spent), #RHS vars (log prices), like Ying's code
              bread=bread_pr_index,
              butter=butter_pr_index,
              cereal=cereal_pr_index,
              chips=chips_pr_index,
              coffee=coffee_pr_index,
              cookies=cookies_pr_index,
              eggs=eggs_pr_index,
              ic=ic_pr_index,
              milk=milk_pr_index,
              oj=oj_pr_index,
              salad=salad_pr_index,
              soda=soda_pr_index,
              soup=soup_pr_index,
              water=water_pr_index,
              yogurt=yogurt_pr_index
              
)

table(df$n)