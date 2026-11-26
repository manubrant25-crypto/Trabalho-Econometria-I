# Econometria I - Trabalho final
# Evasão escolar e sexo, jovens de 15 a 17 anos (PNAD Contínua 2023, 1ª visita)
# José Francisco Lopes Gomes e Manuela Brant Camara

pacotes <- c("PNADcIBGE", "AER", "dplyr", "stargazer", "marginaleffects", "ggplot2", "survey")
faltam <- setdiff(pacotes, rownames(installed.packages()))
if (length(faltam) > 0) install.packages(faltam, repos = "https://cloud.r-project.org")

library(PNADcIBGE)
library(AER)          
library(dplyr)
library(stargazer)
library(marginaleffects)
library(ggplot2)
library(survey)

dir.create("dados/derivados", showWarnings = FALSE, recursive = TRUE)
dir.create("resultados/figuras", showWarnings = FALSE, recursive = TRUE)
dir.create("resultados/tabelas", showWarnings = FALSE, recursive = TRUE)

# Dados ----

arquivo <- "dados/pnadc_2023_visita1.rds"

if (!file.exists(arquivo)) {
  dados <- get_pnadc(year = 2023, interview = 1,
                     vars = c("UF", "V1022", "V2005", "V2007", "V2009", "V2010",
                              "V3002", "VD3004", "VD5008"),
                     design = FALSE, labels = FALSE, deflator = FALSE,
                     savedir = tempdir())
  saveRDS(dados, arquivo)
}

pnad <- readRDS(arquivo) %>%
  mutate(across(c(UF, V1022, V2005, V2007, V2009, V2010, V3002, VD3004, VD5008, V1032),
                as.numeric))

# Amostra ----

# escolaridade do responsável pelo domicílio
resp <- pnad %>%
  filter(V2005 == 1) %>%
  select(UPA, V1008, V1014, educ_resp = VD3004)

stopifnot(!anyDuplicated(resp[c("UPA", "V1008", "V1014")]))

jovens <- pnad %>%
  filter(V2009 >= 15, V2009 <= 17) %>%
  left_join(resp, by = c("UPA", "V1008", "V1014"))

amostra <- jovens %>%
  filter(!is.na(V3002), !is.na(VD3004), !is.na(educ_resp), !is.na(VD5008))

# fluxo da amostra, um critério por vez, por sexo
etapas <- list(
  "Pessoas na base" = pnad,
  "Idade de 15 a 17 anos" = jovens,
  "Frequência e instrução observadas" = filter(jovens, !is.na(V3002), !is.na(VD3004)),
  "Escolaridade do responsável observada" =
    filter(jovens, !is.na(V3002), !is.na(VD3004), !is.na(educ_resp)),
  "Renda observada (amostra final)" = amostra
)
fluxo <- data.frame(
  Etapa = names(etapas),
  N = sapply(etapas, nrow),
  Homens = sapply(etapas, function(d) sum(d$V2007 == 1)),
  Mulheres = sapply(etapas, function(d) sum(d$V2007 == 2)),
  row.names = NULL
)
print(fluxo)
write.csv(fluxo, "resultados/tabelas/fluxo_amostra.csv", row.names = FALSE)

# Variáveis ----

base <- amostra %>%
  mutate(
    # evadiu = não frequenta escola e não terminou o ensino médio
    evasao = as.integer(V3002 == 2 & VD3004 < 5),
    fora_escola = as.integer(V3002 == 2),
    mulher = as.integer(V2007 == 2),
    sexo = factor(mulher, labels = c("Homem", "Mulher")),
    idade = factor(V2009),
    idade_num = V2009,
    raca = case_when(V2010 == 1 ~ "Branca",
                     V2010 %in% c(2, 4) ~ "Preta/Parda",
                     V2010 %in% c(3, 5, 9) ~ "Outras"),   # amarela, indígena, ignorada
    raca = factor(raca, levels = c("Branca", "Preta/Parda", "Outras")),
    regiao = c("Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste")[UF %/% 10],
    regiao = factor(regiao, levels = c("Sudeste", "Norte", "Nordeste", "Sul", "Centro-Oeste")),
    rural = as.integer(V1022 == 2),
    area = ifelse(rural == 1, "Rural", "Urbano"),
    # quem tem renda zero fica com log = 0 e entra na dummy
    renda_zero = as.integer(VD5008 == 0),
    log_renda_pc = ifelse(VD5008 > 0, log(VD5008), 0),
    resp_em = as.integer(educ_resp >= 5),
    resp_f = factor(resp_em, labels = c("sem EM", "com EM")),
    resp_lab = ifelse(resp_em == 1, "Resp. com EM completo", "Resp. sem EM completo"),
    eh_resp = as.integer(V2005 == 1),
    filho_ent = as.integer(V2005 %in% c(4, 5, 6)),
    nivel_ensino = VD3004,
    peso = V1032,
    branca = as.integer(raca == "Branca"),
    preta_parda = as.integer(raca == "Preta/Parda"),
    norte = as.integer(regiao == "Norte"),
    nordeste = as.integer(regiao == "Nordeste"),
    sudeste = as.integer(regiao == "Sudeste"),
    sul = as.integer(regiao == "Sul"),
    centro_oeste = as.integer(regiao == "Centro-Oeste")
  )

stopifnot(!anyNA(base$raca))

# funções ponderadas 
desvio <- function(x, w) sqrt(weighted.mean((x - weighted.mean(x, w))^2, w))

quantil <- function(x, w, p) {
  ordem <- order(x)
  acumulado <- cumsum(w[ordem]) / sum(w)
  sapply(p, function(q) x[ordem][which(acumulado >= q)[1]])
}

correlacao <- function(x, y, w) {
  ok <- complete.cases(x, y)
  cov.wt(cbind(x[ok], y[ok]), wt = w[ok] / sum(w[ok]), cor = TRUE)$cor[1, 2]
}

# erro-padrão de uma proporção, com cluster na UPA
ep_prop <- function(y, w, upa) {
  z <- tapply(w * (y - weighted.mean(y, w)), upa, sum)
  n <- length(z)
  sqrt(n / (n - 1) * sum(z^2)) / sum(w)
}

# quintis de renda calculados entre os próprios jovens
cortes <- quantil(base$VD5008, base$peso, c(0.2, 0.4, 0.6, 0.8))
base$quintil <- factor(findInterval(base$VD5008, cortes, left.open = TRUE) + 1,
                       labels = paste0("Q", 1:5))

base <- base %>%
  select(UPA, peso, evasao, fora_escola, mulher, sexo, idade, idade_num,
         raca, regiao, rural, area, VD5008, renda_zero, log_renda_pc, quintil,
         resp_em, resp_f, resp_lab, eh_resp, filho_ent, nivel_ensino,
         branca, preta_parda, norte, nordeste, sudeste, sul, centro_oeste)

saveRDS(base, "dados/derivados/base_analitica.rds")
write.csv(base %>% select(evasao, fora_escola, mulher, idade_num, raca, regiao, rural,
                          VD5008, renda_zero, log_renda_pc, resp_em),
          "dados/derivados/base_analitica.csv", row.names = FALSE)

# peso dividido pela média: não muda os coeficientes e o glm converge
base$w <- base$peso / mean(base$peso)

# erro-padrão robusto com cluster na UPA; graus de liberdade = nº de UPAs - 1
ep <- function(modelo, dados = base) vcovCL(modelo, cluster = dados$UPA, type = "HC1")
gl <- n_distinct(base$UPA) - 1

# linha do coeficiente de mulher em um MPL
coef_m <- function(modelo, dados = base) {
  coeftest(modelo, vcov. = ep(modelo, dados), df = n_distinct(dados$UPA) - 1)["mulher", ]
}

# Descritivas ----

cat("Evadidos:", sum(base$evasao), "| fora da escola:", sum(base$fora_escola), "\n")
cat("Evasão (%):", 100 * weighted.mean(base$evasao, base$peso),
    "| sem peso:", 100 * mean(base$evasao), "\n")
cat("Fora da escola (%):", 100 * weighted.mean(base$fora_escola, base$peso), "\n")
cat("Outras raças:", sum(base$raca == "Outras"), "| renda zero:", sum(base$renda_zero), "\n")

# nível de instrução de quem evadiu
niveis <- c("Sem instrução", "Fundamental incompleto", "Fundamental completo",
            "Médio incompleto")
comp_evadidos <- base %>%
  filter(evasao == 1) %>%
  group_by(nivel_ensino) %>%
  summarise(n = n(), peso = sum(peso)) %>%
  mutate(nivel = niveis[nivel_ensino],
         pct = 100 * n / sum(n),
         pct_pond = 100 * peso / sum(peso)) %>%
  select(nivel, n, pct, pct_pond)
print(comp_evadidos)
write.csv(comp_evadidos, "resultados/tabelas/composicao_evadidos.csv", row.names = FALSE)

# relação com o responsável
comp_relacao <- base %>%
  mutate(relacao = case_when(eh_resp == 1 ~ "Próprio responsável",
                             filho_ent == 1 ~ "Filho(a) ou enteado(a)",
                             TRUE ~ "Outra relação")) %>%
  group_by(relacao) %>%
  summarise(n = n(),
            pct = 100 * n() / nrow(base),
            pct_pond = 100 * sum(peso) / sum(base$peso),
            evasao_pond = 100 * weighted.mean(evasao, peso))
print(comp_relacao)
write.csv(comp_relacao, "resultados/tabelas/composicao_relacao_resp.csv", row.names = FALSE)

# jovens que são o próprio responsável
proprios <- base %>% filter(eh_resp == 1)
tab_proprios <- table(resp_em = proprios$resp_em, evasao = proprios$evasao)
print(tab_proprios)
write.csv(as.data.frame(tab_proprios), "resultados/tabelas/resp_proprio_resp_em_evasao.csv",
          row.names = FALSE)

# Tabela 1
variaveis <- c("Evasão escolar" = "evasao",
               "Idade (anos)" = "idade_num",
               "Branca" = "branca",
               "Preta/Parda" = "preta_parda",
               "Rural" = "rural",
               "Renda domiciliar pc (R$)" = "VD5008",
               "log(renda pc), se renda > 0" = "log_renda_pc",
               "Renda pc igual a zero" = "renda_zero",
               "Responsável com EM completo" = "resp_em",
               "Norte" = "norte", "Nordeste" = "nordeste", "Sudeste" = "sudeste",
               "Sul" = "sul", "Centro-Oeste" = "centro_oeste")

linha_descritiva <- function(nome) {
  v <- variaveis[[nome]]
  d <- base
  if (v == "log_renda_pc") d <- filter(base, renda_zero == 0)
  homens <- d$mulher == 0
  x <- d[[v]]
  p <- d$peso
  reg <- lm(reformulate("mulher", v), data = d, weights = w)
  teste <- coef_m(reg, d)
  data.frame(Variavel = nome,
             N_obs = nrow(d),
             Media_total = weighted.mean(x, p),
             DP_total = desvio(x, p),
             Media_homem = weighted.mean(x[homens], p[homens]),
             Media_mulher = weighted.mean(x[!homens], p[!homens]),
             Diferenca = teste[["Estimate"]],
             p_valor = teste[["Pr(>|t|)"]])
}

tabela1 <- bind_rows(lapply(names(variaveis), linha_descritiva))
print(tabela1, digits = 3)
write.csv(tabela1, "resultados/tabelas/descritivas.csv", row.names = FALSE)

# correlações
vars_cor <- c("Evasão" = "evasao", "Mulher" = "mulher", "Idade" = "idade_num",
              "Preta/parda" = "preta_parda", "Rural" = "rural",
              "log(renda pc)" = "log_renda_pc", "Resp. com EM" = "resp_em")
d_cor <- base[vars_cor]
d_cor$log_renda_pc[base$renda_zero == 1] <- NA
matriz_cor <- sapply(d_cor, function(x) sapply(d_cor, function(y) correlacao(x, y, base$peso)))
dimnames(matriz_cor) <- list(names(vars_cor), names(vars_cor))
write.csv(round(matriz_cor, 4), "resultados/tabelas/correlacao.csv")

# Gráficos ----

tema <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
cores <- c("Homem" = "#4D4D4D", "Mulher" = "#C0504D")
dodge <- position_dodge(width = 0.8)
quebra <- function(texto, largura) paste(strwrap(texto, largura), collapse = "\n")
legenda <- "Estimativas ponderadas por V1032; IC 95% com erro-padrão por conglomerados (UPA)."

salvar <- function(grafico, nome, largura = 7, altura = 4.2) {
  ggsave(paste0("resultados/figuras/", nome, ".png"), grafico,
         width = largura, height = altura, dpi = 200)
}

taxa_por <- function(...) {
  base %>%
    group_by(...) %>%
    summarise(n = n(),
              p = weighted.mean(evasao, peso),
              se = ep_prop(evasao, peso, UPA),
              .groups = "drop") %>%
    mutate(evasao = 100 * p,
           ci_l = 100 * pmax(0, p - 1.96 * se),
           ci_u = 100 * pmin(1, p + 1.96 * se))
}

g1 <- ggplot(taxa_por(sexo), aes(sexo, evasao, fill = sexo)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), width = 0.15) +
  geom_text(aes(y = ci_u, label = sprintf("%.1f%%", evasao)), vjust = -0.6, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = cores, guide = "none") +
  labs(x = NULL, y = "Taxa de evasão (%)", title = "Evasão escolar por sexo",
       subtitle = "Jovens de 15 a 17 anos, PNADC 2023", caption = quebra(legenda, 60)) +
  tema
salvar(g1, "fig1_evasao_sexo", 4.5, 4)

barras <- function(dados, x) {
  ggplot(dados, aes({{ x }}, evasao, fill = sexo)) +
    geom_col(position = dodge, width = 0.7) +
    geom_errorbar(aes(ymin = ci_l, ymax = ci_u), position = dodge, width = 0.2) +
    scale_fill_manual(values = cores, name = NULL)
}


salvar_painel <- function(grafico, nome) {
  ggsave(paste0("resultados/figuras/", nome, ".png"), grafico,
         width = 3.4, height = 2.04, dpi = 300)
}
tema_painel <- tema +
  theme(legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        legend.key.size = unit(8, "pt"), plot.margin = margin(2, 6, 2, 2))

g2 <- barras(taxa_por(idade, sexo), idade) +
  labs(x = "Idade", y = "Evasão (%)") +
  tema_painel
salvar_painel(g2, "fig2_evasao_idade_sexo")

g3 <- barras(taxa_por(regiao, sexo), regiao) +
  scale_x_discrete(labels = c("Centro-Oeste" = "Centro-\nOeste")) +
  labs(x = NULL, y = "Evasão (%)") +
  tema_painel
salvar_painel(g3, "fig3_evasao_regiao_sexo")

g4 <- ggplot(taxa_por(quintil, sexo), aes(quintil, evasao, color = sexo, group = sexo)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), width = 0.12) +
  scale_color_manual(values = cores, name = NULL) +
  labs(x = "Quintil de renda (Q1 = mais pobre)", y = "Evasão (%)") +
  tema_painel
salvar_painel(g4, "fig4_evasao_renda_sexo")

g5 <- barras(taxa_por(resp_lab, area, sexo), area) +
  facet_wrap(~resp_lab) +
  labs(x = NULL, y = "Taxa de evasão (%)",
       title = "Evasão por área, escolaridade do responsável e sexo",
       caption = quebra(legenda, 95)) +
  tema
salvar(g5, "fig5_evasao_area_resp_sexo")

com_renda <- base %>%
  filter(renda_zero == 0) %>%
  mutate(situacao = factor(evasao, labels = c("Não evadiu", "Evadiu")))

g6 <- ggplot(com_renda, aes(log_renda_pc, weight = peso / sum(peso), fill = situacao)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#7F7F7F", "#C0504D"), name = NULL) +
  labs(x = "log(renda domiciliar per capita)", y = "Densidade",
       title = "Distribuição da renda por situação de evasão",
       subtitle = "Exclui renda zero. Densidades ponderadas por V1032") +
  tema
salvar(g6, "fig6_densidade_renda")

# renda em reais e em log 
p_acima <- 100 * weighted.mean(com_renda$VD5008 > 6000, com_renda$peso)
cat("Renda acima de R$ 6.000 (%):", p_acima, "\n")

paineis <- c("Renda pc (R$)", "log(renda pc)")
dois <- bind_rows(
  com_renda %>% filter(VD5008 <= 6000) %>% mutate(x = VD5008, painel = paineis[1]),
  com_renda %>% mutate(x = log_renda_pc, painel = paineis[2])
) %>%
  mutate(painel = factor(painel, levels = paineis))

g9 <- ggplot(dois, aes(x, weight = peso / sum(peso), fill = situacao)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~painel, scales = "free") +
  scale_x_continuous(breaks = function(l) if (l[2] > 100) c(0, 3000, 6000) else c(2, 4, 6, 8, 10)) +
  scale_fill_manual(values = c("#7F7F7F", "#C0504D"), name = NULL) +
  labs(x = NULL, y = "Densidade") +
  tema_painel +
  theme(axis.text.y = element_blank())
salvar_painel(g9, "fig9_assimetria_renda")

# Modelos ----

# trabalho, gravidez e estado civil ficam de fora por serem mediadores
principal <- evasao ~ mulher + idade + raca + regiao + rural + resp_em + log_renda_pc + renda_zero

# confere multicolinearidade perfeita
X <- model.matrix(principal, data = base)
cat("Posto:", qr(X)$rank, "de", ncol(X), "colunas\n")

# MPL
m1 <- lm(evasao ~ mulher, data = base, weights = w)
m2 <- lm(evasao ~ mulher + idade + raca + regiao + rural, data = base, weights = w)
m3 <- lm(evasao ~ mulher + idade + raca + regiao + rural + resp_em, data = base, weights = w)
m4 <- lm(principal, data = base, weights = w)

# multicolinearidade imperfeita (GVIF para os fatores)
print(vif(m4))

# Probit e Logit (o aviso de "non-integer #successes" vem dos pesos)
# critério de convergência mais apertado que o padrão do glm
conv <- glm.control(epsilon = 1e-12, maxit = 100)
probit <- glm(principal, family = binomial(link = "probit"), data = base, weights = w, control = conv)
logit <- glm(principal, family = binomial(link = "logit"), data = base, weights = w, control = conv)

# Tabela 3
testes_coef <- list(coeftest(m1, vcov. = ep(m1), df = gl),
                    coeftest(m2, vcov. = ep(m2), df = gl),
                    coeftest(m3, vcov. = ep(m3), df = gl),
                    coeftest(m4, vcov. = ep(m4), df = gl),
                    coeftest(probit, vcov. = ep(probit)),
                    coeftest(logit, vcov. = ep(logit)))

stargazer(m1, m2, m3, m4, probit, logit,
          se = lapply(testes_coef, function(t) t[, 2]),
          p = lapply(testes_coef, function(t) t[, 4]),
          type = "text", digits = 4, star.cutoffs = c(0.1, 0.05, 0.01),
          keep.stat = c("n", "rsq"), out = "resultados/tabelas/tabela3.txt")

# efeito marginal médio
ame <- function(modelo, grupo = TRUE) {
  avg_slopes(modelo, variables = "mulher", by = grupo, wts = "peso", vcov = ep(modelo))
}

ame_probit <- ame(probit)
ame_logit <- ame(logit)
print(ame_probit)
print(ame_logit)

mpl <- list(m1, m2, m3, m4)
est <- t(sapply(mpl, coef_m))
ic <- t(sapply(mpl, function(m) coefci(m, "mulher", vcov. = ep(m), df = gl)))

coef_mulher <- data.frame(
  Modelo = c("MPL sem controles", "MPL + demográficos/região", "MPL + escolaridade do resp.",
             "MPL principal (+ renda)", "Probit (AME)", "Logit (AME)"),
  Coef = c(est[, 1], ame_probit$estimate, ame_logit$estimate),
  EP = c(est[, 2], ame_probit$std.error, ame_logit$std.error),
  IC_inf = c(ic[, 1], ame_probit$conf.low, ame_logit$conf.low),
  IC_sup = c(ic[, 2], ame_probit$conf.high, ame_logit$conf.high)
)
print(coef_mulher, digits = 3)
write.csv(coef_mulher, "resultados/tabelas/coef_mulher_modelos.csv", row.names = FALSE)

# Testes e diagnósticos ----

# as especificações com interação partem sempre da principal
f_regiao <- update(principal, . ~ . + mulher:regiao)
mpl_regiao <- lm(f_regiao, data = base, weights = w)
probit_regiao <- glm(f_regiao, family = binomial(link = "probit"), data = base, weights = w, control = conv)

# teste F conjunto com a variância por cluster
teste_f <- function(modelo, termos) {
  linearHypothesis(modelo, matchCoefs(modelo, termos), vcov. = ep(modelo),
                   test = "F", error.df = gl)
}

linha_teste <- function(modelo, termos, hipotese, tipo = "MPL") {
  t <- teste_f(modelo, termos)
  data.frame(Hipotese = hipotese, Modelo = tipo, Estatistica = t$F[2],
             p_valor = t$`Pr(>F)`[2], gl1 = t$Df[2], gl2 = gl)
}

testes <- rbind(
  linha_teste(m4, "^idade", "Idade (todas dummies = 0)"),
  linha_teste(m4, "^raca", "Raça (todas dummies = 0)"),
  linha_teste(m4, "^regiao", "Região (todas dummies = 0)"),
  linha_teste(m4, "renda", "Renda (log + zero = 0)"),
  linha_teste(mpl_regiao, "mulher:regiao", "Interações mulher x região = 0"),
  linha_teste(probit_regiao, "mulher:regiao", "Interações mulher x região = 0", "Probit")
)
print(testes, digits = 3, row.names = FALSE)
write.csv(testes, "resultados/tabelas/testes_conjuntos.csv", row.names = FALSE)

cat("Número de UPAs:", n_distinct(base$UPA), "\n")

# previsões do MPL fora de [0,1]
previsto <- predict(m4)
cat("Mínimo:", min(previsto), "| máximo:", max(previsto),
    "| % negativas:", 100 * mean(previsto < 0), "\n")

# matriz de confusão do probit (sem pesos)
prob <- predict(probit, type = "response")
print(table(observado = base$evasao, previsto = as.integer(prob > mean(base$evasao))))
cat("Acerto com corte 0,5:", mean((prob > 0.5) == base$evasao), "\n")

# calibração por decil
calib <- base %>%
  mutate(prob = prob,
         decil = cut(prob, quantile(prob, seq(0, 1, 0.1)), include.lowest = TRUE)) %>%
  group_by(decil) %>%
  summarise(previsto = 100 * weighted.mean(prob, peso), observado = 100 * weighted.mean(evasao, peso))

g7 <- ggplot(calib, aes(previsto, observado)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 2.5, color = "#C0504D") +
  labs(x = "Probabilidade prevista média (%)", y = "Taxa observada (%)",
       title = "Calibração do Probit por decil de risco",
       caption = quebra("Médias ponderadas por V1032; decis da probabilidade prevista (amostrais).", 70)) +
  tema
salvar(g7, "fig7_calibracao_probit", 5, 4.2)

# teste de ligação 
base$xb_probit <- predict(probit, type = "link")
base$xb_logit <- predict(logit, type = "link")
base$ajustado <- previsto

liga_probit <- glm(update(principal, . ~ . + I(xb_probit^2)),
                   family = binomial(link = "probit"), data = base, weights = w, control = conv)
liga_logit <- glm(update(principal, . ~ . + I(xb_logit^2)),
                  family = binomial(link = "logit"), data = base, weights = w, control = conv)
print(coeftest(liga_probit, vcov. = ep(liga_probit)))
print(coeftest(liga_logit, vcov. = ep(liga_logit)))

# RESET
reset <- lm(update(principal, . ~ . + I(ajustado^2) + I(ajustado^3)), data = base, weights = w)
print(teste_f(reset, "ajustado"))

# erro-padrão com os pesos replicados do IBGE
replicas <- amostra %>%
  select(matches("^V1032[0-9]{3}$")) %>%
  mutate(across(everything(), as.numeric))
# 200 réplicas, na mesma ordem de linhas da base
stopifnot(ncol(replicas) == 200, !anyNA(replicas), identical(amostra$V1032, base$peso))
# bootstrap com mse = TRUE, a mesma configuração do pacote PNADcIBGE
desenho <- svrepdesign(data = base, repweights = as.matrix(replicas), weights = ~peso,
                       type = "bootstrap", mse = TRUE)
mpl_replicas <- svyglm(principal, design = desenho)
cat("EP por UPA:", coef_m(m4)[["Std. Error"]],
    "| EP com pesos replicados:", summary(mpl_replicas)$coefficients["mulher", "Std. Error"], "\n")

# Robustez ----

d17 <- filter(base, idade_num == 17)
d1516 <- filter(base, idade_num < 17)
d_filhos <- filter(base, filho_ent == 1)
d_semresp <- filter(base, eh_resp == 0)
d_semzero <- filter(base, renda_zero == 0)
d_seminstr <- filter(base, nivel_ensino != 1)

fora <- lm(update(principal, fora_escola ~ .), data = base, weights = w)
# com 17 anos só, a idade é constante
so17 <- lm(update(principal, . ~ . - idade), data = d17, weights = w)
de15a16 <- lm(principal, data = d1516, weights = w)
filhos <- lm(principal, data = d_filhos, weights = w)
sem_pesos <- lm(principal, data = base)
sem_resp <- lm(principal, data = d_semresp, weights = w)
sem_zero <- lm(update(principal, . ~ . - renda_zero), data = d_semzero, weights = w)
sem_instrucao <- lm(principal, data = d_seminstr, weights = w)

# Tabela 4 
modelos <- list("Principal" = list(m4, base),
                "Sem renda" = list(m3, base),
                "Fora da escola" = list(fora, base),
                "17 anos" = list(so17, d17),
                "15-16 anos" = list(de15a16, d1516),
                "Filho/enteado" = list(filhos, d_filhos),
                "Sem pesos" = list(sem_pesos, base),
                "Sem o resp." = list(sem_resp, d_semresp),
                "Sem renda zero" = list(sem_zero, d_semzero),
                "Sem instrução" = list(sem_instrucao, d_seminstr))

tabela4 <- bind_rows(lapply(names(modelos), function(nome) {
  m <- modelos[[nome]][[1]]
  t <- coef_m(m, modelos[[nome]][[2]])
  data.frame(Especificacao = nome, Coef = t[[1]], EP = t[[2]], p_valor = t[[4]],
             N = nobs(m), R2 = summary(m)$r.squared)
}))
print(tabela4, digits = 3)
write.csv(tabela4, "resultados/tabelas/robustez.csv", row.names = FALSE)

cat("Excluídos por 'sem instrução':", sum(base$nivel_ensino == 1), "\n")

# Heterogeneidade ----

# por região
het_regiao <- as.data.frame(ame(probit_regiao, "regiao")) %>%
  select(regiao, estimate, std.error, p.value, conf.low, conf.high)
print(het_regiao, digits = 3)
write.csv(het_regiao, "resultados/tabelas/ame_mulher_por_regiao.csv", row.names = FALSE)

pontos <- bind_rows(
  coef_mulher %>%
    transmute(grupo = "Modelos", rotulo = Modelo,
              b = 100 * Coef, lo = 100 * IC_inf, hi = 100 * IC_sup),
  het_regiao %>%
    transmute(grupo = "AME por região (Probit)", rotulo = as.character(regiao),
              b = 100 * estimate, lo = 100 * conf.low, hi = 100 * conf.high)
)
pontos$rotulo <- factor(pontos$rotulo, levels = rev(unique(pontos$rotulo)))

g8 <- ggplot(pontos, aes(b, rotulo)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), color = "#C0504D") +
  facet_grid(grupo ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Diferença Mulher - Homem na probabilidade de evasão (p.p.), IC 95%", y = NULL,
       title = "Associação entre ser mulher e evasão",
       caption = quebra("Estimativas ponderadas por V1032; EP por conglomerados (UPA).", 95)) +
  tema
salvar(g8, "fig8_coef_mulher", 7, 5)

# por idade e por escolaridade do responsável 
f_idade <- update(principal, . ~ . + mulher:idade)
f_resp <- update(principal, . ~ . + mulher:resp_em)

probit_idade <- glm(f_idade, family = binomial(link = "probit"), data = base, weights = w, control = conv)
probit_resp <- glm(update(principal, . ~ . - resp_em + resp_f + mulher:resp_f),
                   family = binomial(link = "probit"), data = base, weights = w, control = conv)
print(ame(probit_idade, "idade"))
print(ame(probit_resp, "resp_f"))

# interações no MPL
mpl_idade <- lm(f_idade, data = base, weights = w)
mpl_resp <- lm(f_resp, data = base, weights = w)
print(teste_f(mpl_idade, "mulher:idade"))
print(teste_f(mpl_resp, "mulher:resp_em"))

# Poisson, para testar se a razão mulher/homem muda entre os grupos
pois_idade <- glm(f_idade, family = poisson, data = base, weights = w)
pois_resp <- glm(f_resp, family = poisson, data = base, weights = w)
pois_regiao <- glm(f_regiao, family = poisson, data = base, weights = w)

print(teste_f(pois_idade, "mulher:idade"))
print(teste_f(pois_resp, "mulher:resp_em"))
print(teste_f(pois_regiao, "mulher:regiao"))

b <- coef(pois_idade)
cat("Razões por idade:", exp(b["mulher"] + c(0, b["mulher:idade16"], b["mulher:idade17"])), "\n")
b <- coef(pois_resp)
cat("Razões por responsável:", exp(b["mulher"] + c(0, b["mulher:resp_em"])), "\n")

# teste de Wald para igualdade dos AME do probit entre os grupos
igualdade_ame <- function(modelo, grupo) {
  a <- ame(modelo, grupo)
  k <- nrow(a)
  R <- cbind(1, -diag(k - 1))
  d <- R %*% a$estimate
  estat <- as.numeric(t(d) %*% solve(R %*% vcov(a) %*% t(R)) %*% d) / (k - 1)
  p <- pf(estat, k - 1, gl, lower.tail = FALSE)
  cat(grupo, ": F =", estat, " p =", p, "\n")
}

igualdade_ame(probit_idade, "idade")
igualdade_ame(probit_resp, "resp_f")
igualdade_ame(probit_regiao, "regiao")

# todos os glm precisam ter convergido
modelos_glm <- list(probit, logit, probit_regiao, probit_idade, probit_resp,
                    liga_probit, liga_logit, pois_idade, pois_resp, pois_regiao)
convergiu <- sapply(modelos_glm, function(m) m$converged)
cat("Modelos glm que convergiram:", sum(convergiu), "de", length(convergiu), "\n")
stopifnot(all(convergiu))

# taxas de evasão citadas no texto
taxa <- function(...) {
  base %>%
    group_by(...) %>%
    summarise(n = n(), evasao = 100 * weighted.mean(evasao, peso), .groups = "drop")
}
print(taxa(sexo))
print(taxa(resp_f, sexo))
print(taxa(idade_num, sexo))
print(taxa(regiao, sexo))
print(taxa(quintil, sexo))
print(taxa(filho_ent, sexo))

writeLines(capture.output(sessionInfo()), "resultados/sessionInfo.txt")
