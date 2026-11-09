# compositeSS graphical interface
# Launch with compositeSS::run_calculator() (or compositeSS::cs_app()).
library(shiny)
library(compositeSS)

num <- function(id, label, value, min = NA, max = NA, step = NA)
  numericInput(id, label, value, min = min, max = max, step = step)

safe <- function(expr) tryCatch(expr, error = function(e) {
  validate(need(FALSE, conditionMessage(e)))
})

# ---------------------------------------------------------------- UI ----
design_inputs <- tagList(
  h4("Estimand"),
  radioButtons("estimand", NULL,
               c("Shih-Quan composite test" = "cs",
                 "Expected benefit per randomized patient (theta)" = "theta")),
  h4("Planning values"),
  fluidRow(column(6, num("pA", "Completion, A", 0.65, 0.01, 0.99, 0.01)),
           column(6, num("pB", "Completion, B", 0.50, 0.01, 0.99, 0.01))),
  conditionalPanel("input.efficacy != 'binary' || input.estimand == 'theta'",
    fluidRow(column(6, num("muA", "Mean improvement, completers A", 21)),
             column(6, num("muB", "Mean improvement, completers B", 15))),
    num("sigma", "SD among completers", 18, 0.1)),
  conditionalPanel("input.estimand == 'cs'",
    selectInput("efficacy", "Completers' efficacy outcome",
                c("Continuous" = "continuous", "Response (binary)" = "binary")),
    conditionalPanel("input.efficacy == 'binary'",
      fluidRow(column(6, num("respA", "Response, A", 0.55, 0.01, 0.99, 0.01)),
               column(6, num("respB", "Response, B", 0.40, 0.01, 0.99, 0.01)))),
    h4("Test"),
    selectInput("method", "Combination",
                c("Weighted Z" = "wz", "Fisher (equal weights)" = "fisher")),
    conditionalPanel("input.method == 'wz' && input.efficacy == 'binary'",
      checkboxInput("aligned", "Estimand-aligned weights (test the probability of completing and responding)")),
    conditionalPanel("input.method == 'wz' && !(input.aligned && input.efficacy == 'binary')",
      sliderInput("wd", "Weight on completion, w(d)  [w(e) = 1 - w(d)]",
                  0, 1, 0.5, 0.05)),
    checkboxInput("contingent", "Contingent test (stop if A worse on completion)"),
    conditionalPanel("input.contingent",
      fluidRow(column(6, num("gate", "Gate c", 0, step = 0.1)),
               column(6, selectInput("cal", "Critical value",
                                     c("recalibrated", "usual"))))),
    selectInput("gamma", "Covariance gamma under H1",
                c("Computed (delta method)" = "compute", "Zero" = "zero")),
    conditionalPanel("input.method == 'wz'",
      fluidRow(column(6, num("margin_d", "NI margin, completion", 0, 0, 0.5, 0.01)),
               column(6, num("margin_e", "NI margin, efficacy", 0, 0))))),
  conditionalPanel("input.estimand == 'theta'",
    fluidRow(column(6, num("s", "Dropout score s", 0)),
             column(6, num("margin_t", "NI margin (clinical units)", 0, 0)))),
  h4("Error rates and allocation"),
  fluidRow(column(4, num("alpha", "One-sided alpha", 0.025, 0.001, 0.25, 0.005)),
           column(4, num("power", "Target power", 0.80, 0.5, 0.99, 0.01)),
           column(4, num("r", "nB / nA", 1, 0.1, 10, 0.1)))
)

ui <- navbarPage(
  "compositeSS",
  tabPanel("Sample size",
    sidebarLayout(
      sidebarPanel(design_inputs, width = 4),
      mainPanel(
        h4("Required sample size"), verbatimTextOutput("ss_txt"),
        plotOutput("curve", height = 320),
        h4("Power at the required size"), tableOutput("pw_tab"),
        h4("Estimand (ICH E9(R1))"), verbatimTextOutput("estimand_txt"),
        width = 8))),

  tabPanel("Compare designs",
    sidebarLayout(
      sidebarPanel(
        p("Sample sizes per arm for the designs compared in Table 6, using the",
          "planning values, alpha and power on the Sample size tab (1:1 allocation)."),
        selectInput("cmp_gamma", "gamma", c("Zero (as in Table 6)" = "zero",
                                            "Computed" = "compute")),
        actionButton("cmp_t6", "Load Table 6 scenarios"),
        width = 3),
      mainPanel(tableOutput("cmp_tab"), plotOutput("cmp_plot", height = 360),
                h4("Table 6 of the proposal"), tableOutput("t6_tab"), width = 9))),

  tabPanel("Sensitivity & assurance",
    sidebarLayout(
      sidebarPanel(
        p("Uses the design on the Sample size tab."),
        num("sens_n", "Patients on A (blank = required size)", NA, 2),
        h4("Grid of true values"),
        sliderInput("sens_pA", "True completion A", 0.2, 0.95, c(0.50, 0.80), 0.01),
        sliderInput("sens_pB", "True completion B", 0.2, 0.95, c(0.35, 0.65), 0.01),
        h4("Assurance priors"),
        p("Beta priors given by a mean and an effective sample size; normal",
          "prior for the completers' difference."),
        fluidRow(column(6, num("pr_pA_m", "Mean, completion A", 0.65, 0.01, 0.99)),
                 column(6, num("pr_pA_n", "Eff. n", 100, 1))),
        fluidRow(column(6, num("pr_pB_m", "Mean, completion B", 0.50, 0.01, 0.99)),
                 column(6, num("pr_pB_n", "Eff. n", 100, 1))),
        fluidRow(column(6, num("pr_d_m", "Mean, difference", 6)),
                 column(6, num("pr_d_sd", "SD, difference", 2, 0))),
        num("pr_nsim", "Prior draws", 2000, 100),
        actionButton("go_assur", "Compute assurance"),
        width = 4),
      mainPanel(plotOutput("heat", height = 380), plotOutput("delta_curve", height = 280),
                h4("Assurance"), verbatimTextOutput("assur_txt"),
                plotOutput("assur_hist", height = 250), width = 8))),

  tabPanel("Simulation",
    sidebarLayout(
      sidebarPanel(
        p("Simulates the composite test on the Sample size tab",
          "(continuous efficacy, superiority)."),
        fluidRow(column(6, num("sim_n", "Patients on A (blank = required)", NA, 2)),
                 column(6, num("sim_nsim", "Simulated trials", 500, 10))),
        checkboxInput("sim_null", "Simulate under the null hypothesis"),
        radioButtons("sim_dropout", "Dropout model",
                     c("Unrelated to outcome (Appendix E)" = "random",
                       "Appendix C longitudinal model" = "appendixC")),
        conditionalPanel("input.sim_dropout == 'appendixC'",
          h4("Part 1: trajectories"),
          fluidRow(column(6, num("k", "Visits k", 6, 2, 12, 1)),
                   column(6, selectInput("shape", "Mean trajectory",
                                         c("quadratic", "linear", "flat")))),
          selectInput("covtype", "Covariance",
                      c("Compound symmetry" = "cs",
                        "Random intercept and slope" = "ris")),
          conditionalPanel("input.covtype == 'cs'",
            fluidRow(column(6, num("csd", "SD per visit", 18, 0.1)),
                     column(6, num("rho", "Correlation rho", 0.5, 0, 0.95, 0.05)))),
          conditionalPanel("input.covtype == 'ris'",
            fluidRow(column(6, num("sd_int", "SD intercept", 12, 0)),
                     column(6, num("sd_slope", "SD slope", 10, 0))),
            fluidRow(column(6, num("cor_is", "Corr(int, slope)", 0, -0.95, 0.95, 0.05)),
                     column(6, num("sd_resid", "Residual SD", 10, 0.1)))),
          h4("Part 2: dropout"),
          num("lambda", "Outcome-unrelated dropout per interval", 0.02, 0, 0.5, 0.01),
          fluidRow(column(6, num("betaA", "beta, A", -0.08, step = 0.01)),
                   column(6, num("betaB", "beta, B", -0.08, step = 0.01))),
          helpText("beta < 0: less improvement, more dropout (lack of response)."),
          selectInput("depends", "Dropout depends on",
                      c("Current (unobserved) value: MNAR" = "current",
                        "Last observed value: MAR" = "previous")),
          num("n_cal", "Calibration sample", 20000, 2000),
          actionButton("go_preview", "Calibrate and preview model"),
          h4("Comparators"),
          checkboxGroupInput("comparators", NULL,
                             c("MMRM" = "mmrm",
                               "Placebo-based multiple imputation" = "pmi")),
          fluidRow(column(6, num("M", "Imputations", 20, 2)),
                   column(6, selectInput("engine", "MMRM engine",
                                         c("nlme", "mmrm"))))),
        num("sim_seed", "Seed", 2026, 1),
        actionButton("go_sim", "Run simulation", class = "btn-primary"),
        width = 4),
      mainPanel(plotOutput("appc_plot", height = 330),
                verbatimTextOutput("appc_txt"),
                h4("Simulation results"), verbatimTextOutput("sim_txt"),
                tableOutput("cmp_sim_tab"), width = 8))),

  tabPanel("Analyze data",
    sidebarLayout(
      sidebarPanel(
        fileInput("file", "Trial data (CSV, one row per randomized patient)",
                  accept = ".csv"),
        uiOutput("col_ui"),
        h4("Analysis"),
        selectInput("an_method", "Combination", c("Weighted Z" = "wz",
                                                  "Fisher" = "fisher")),
        conditionalPanel("input.an_method == 'wz'",
          checkboxInput("an_aligned", "Estimand-aligned weights (computed from the data)")),
        conditionalPanel("input.an_method == 'wz' && !input.an_aligned",
          sliderInput("an_wd", "Weight on completion", 0, 1, 0.5, 0.05)),
        checkboxInput("an_cont", "Contingent test"),
        conditionalPanel("input.an_cont",
          num("an_c", "Gate c", 0, step = 0.1)),
        fluidRow(column(6, num("an_alpha", "One-sided alpha", 0.025, 0.001, 0.25)),
                 column(6, num("an_s", "Dropout score s", 0))),
        selectInput("an_psi", "psi", c("information", "clinical")),
        num("an_boot", "Bootstrap resamples (0 = none)", 0, 0),
        actionButton("go_an", "Analyze", class = "btn-primary"),
        width = 4),
      mainPanel(verbatimTextOutput("an_txt"), tableOutput("an_head"), width = 8))),

  tabPanel("About",
    fluidPage(
      h3("compositeSS"),
      p("Design and sample size methods for trials in which dropout is part of",
        "the outcome (composite strategy of ICH E9(R1)): the expected benefit",
        "per randomized patient and the Shih-Quan composite test with unequal",
        "weights, the standardized composite effect size psi and contingent",
        "testing."),
      p("Simulation uses either dropout unrelated to outcome (Appendix E) or",
        "the three-part longitudinal model of Appendix C, with MMRM and",
        "placebo-based multiple imputation as comparators judged against the",
        "treatment policy effect."),
      p("Improvement is coded so that positive values favour the patient.")))
)

# ------------------------------------------------------------ server ----
server <- function(input, output, session) {

  cs_des <- reactive({
    wd <- if (input$method == "fisher") 0.5 else input$wd
    binary <- input$efficacy == "binary"
    aligned <- binary && input$method == "wz" && isTRUE(input$aligned)
    safe(cs_design(input$pA, input$pB,
                   delta = if (binary) NULL else input$muA - input$muB,
                   sigma = if (binary) NULL else input$sigma,
                   r = input$r,
                   weights = if (aligned) "aligned" else c(wd, 1 - wd),
                   method = input$method,
                   contingent = input$contingent,
                   c = if (is.null(input$gate)) 0 else input$gate,
                   calibration = input$cal, alpha = input$alpha,
                   gamma = input$gamma, efficacy = input$efficacy,
                   respA = input$respA, respB = input$respB,
                   margin_d = if (input$method == "wz") input$margin_d else 0,
                   margin_e = if (input$method == "wz") input$margin_e else 0))
  })
  th_des <- reactive(safe(theta_design(input$pA, input$pB, input$muA, input$muB,
                                       input$sigma, s = input$s, r = input$r,
                                       alpha = input$alpha, margin = input$margin_t)))
  des <- reactive(if (input$estimand == "theta") th_des() else cs_des())
  ss <- reactive({
    d <- des()
    safe(if (inherits(d, "theta_design")) theta_sample_size(d, input$power)
         else cs_sample_size(d, input$power))
  })

  # ---- Sample size tab
  output$ss_txt <- renderPrint({
    x <- ss(); d <- des()
    cat(sprintf("nA = %d, nB = %d, total = %d   (power %.3f)\n",
                x$nA, x$nB, x$total, x$power))
    if (inherits(d, "theta_design")) {
      dc <- theta_decompose(d)
      cat(sprintf("theta = %.3f = completion part %.3f + completers part %.3f\n",
                  dc$theta, dc$completion, dc$completers))
    } else {
      if (isTRUE(d$aligned)) {
        cat(sprintf("Estimand-aligned weights: %.3f / %.3f (the test targets the probability of completing and responding)\n",
                    d$weights[1], d$weights[2]))
      } else if (d$method == "wz") {
        ow <- cs_opt_weights(d, x$nA)
        cat(sprintf("Power-maximizing weights: %.2f / %.2f  (power %.3f vs %.3f with chosen weights)\n",
                    ow$weights[1], ow$weights[2], ow$power_optimal, ow$power_design))
      }
      cat(sprintf("Critical value: %.4f (%s scale)\n", x$critical$k, x$critical$scale))
      cat(sprintf("gamma = %.4f\n", cs_components(d, x$nA)$gamma))
      if (d$efficacy == "continuous")
        cat(sprintf("psi = %.3f (information weights), %.3f (clinical weights)\n",
                    cs_psi(d)$psi, cs_psi(d, "clinical")$psi))
    }
    cat("No inflation for dropout is needed: dropout is part of the endpoint.\n")
  })
  output$curve <- renderPlot({
    n_hi <- max(50, ceiling(2 * ss()$nA))
    cs_power_curve(list(Design = des()), n = unique(round(seq(5, n_hi, length.out = 80))),
                   target = input$power, main = "Power by patients per arm")
    abline(v = ss()$nA, lty = 3)
  })
  output$pw_tab <- renderTable({
    d <- des(); n <- ss()$nA
    if (inherits(d, "theta_design"))
      data.frame(nA = n, power = theta_power(d, n))
    else as.data.frame(cs_power(d, n))
  }, digits = 4)
  output$estimand_txt <- renderPrint(invisible(cs_estimand_statement(des())))

  # ---- Compare designs tab
  observeEvent(input$cmp_t6, {
    updateNumericInput(session, "pA", value = 0.65)
    updateNumericInput(session, "pB", value = 0.50)
    updateNumericInput(session, "muA", value = 21)
    updateNumericInput(session, "muB", value = 15)
    updateNumericInput(session, "sigma", value = 18)
    updateNumericInput(session, "alpha", value = 0.025)
    updateNumericInput(session, "power", value = 0.80)
  })
  cmp <- reactive(safe(cs_compare_designs(input$pA, input$pB, input$muA - input$muB,
                                          input$sigma, muB = input$muB, s = input$s,
                                          alpha = input$alpha, power = input$power,
                                          gamma = input$cmp_gamma)))
  output$cmp_tab <- renderTable(transform(cmp(), n_per_arm = as.integer(n_per_arm)), digits = 3)
  output$cmp_plot <- renderPlot({
    x <- cmp()
    op <- par(mar = c(4, 22, 1, 1)); on.exit(par(op))
    barplot(rev(x$n_per_arm), names.arg = rev(x$design), horiz = TRUE, las = 1,
            xlab = "Patients per arm", cex.names = 0.8, col = "grey70")
  })
  output$t6_tab <- renderTable(cs_table6(), digits = 0)

  # ---- Sensitivity tab
  sens_n <- reactive(if (is.na(input$sens_n)) ss()$nA else input$sens_n)
  output$heat <- renderPlot({
    d <- des()
    pa <- seq(input$sens_pA[1], input$sens_pA[2], length.out = 21)
    pb <- seq(input$sens_pB[1], input$sens_pB[2], length.out = 21)
    n <- sens_n()
    z <- if (inherits(d, "theta_design")) {
      outer(pa, pb, Vectorize(function(a, b) {
        d2 <- d; d2$pA <- a; d2$pB <- b; theta_power(d2, n) }))
    } else {
      g <- cs_power_grid(d, n, pA = pa, pB = pb)
      matrix(g$power, length(pa), length(pb))
    }
    image(pa, pb, z, zlim = c(0, 1), col = hcl.colors(40, "Blues 3", rev = TRUE),
          xlab = "True completion, A", ylab = "True completion, B",
          main = sprintf("Power at nA = %d", n))
    contour(pa, pb, z, levels = c(0.5, 0.7, 0.8, 0.9), add = TRUE)
    points(d$pA, d$pB, pch = 19)
  })
  output$delta_curve <- renderPlot({
    d <- des(); n <- sens_n()
    eff0 <- input$muA - input$muB
    eff <- seq(0, 2 * max(abs(eff0), 1), length.out = 50)
    pw <- vapply(eff, function(e) {
      if (inherits(d, "theta_design")) { d2 <- d; d2$muA <- d$muB + e; theta_power(d2, n) }
      else if (d$efficacy == "binary") NA_real_
      else cs_power_grid(d, n, delta = e)$power
    }, numeric(1))
    plot(eff, pw, type = "l", lwd = 2, ylim = c(0, 1),
         xlab = "True completers' difference", ylab = "Power",
         main = "Power if the completers' effect differs from the planning value")
    abline(v = eff0, lty = 3); abline(h = input$power, col = "grey50")
  })
  assur <- eventReactive(input$go_assur, {
    d <- des()
    bp <- function(m, n) c(m * n, (1 - m) * n)
    safe(cs_assurance(d, sens_n(), prior_pA = bp(input$pr_pA_m, input$pr_pA_n),
                      prior_pB = bp(input$pr_pB_m, input$pr_pB_n),
                      prior_delta = if (input$pr_d_sd > 0)
                        c(input$pr_d_m, input$pr_d_sd) else NULL,
                      nsim = input$pr_nsim, seed = 1))
  })
  output$assur_txt <- renderPrint({
    a <- assur()
    cat(sprintf("Assurance at nA = %d: %.3f   (power at planning values %.3f)\n",
                sens_n(), a$assurance, a$power_at_planning_values))
  })
  output$assur_hist <- renderPlot({
    a <- assur()
    hist(a$draws$power, breaks = 30, main = "Power over prior draws",
         xlab = "Power", col = "grey80")
    abline(v = a$assurance, lwd = 2, col = 2)
  })

  # ---- Simulation tab
  appc_model <- reactive(safe(cs_appc_model(
    k = input$k, shape = input$shape, cov = input$covtype, sd = input$csd,
    rho = input$rho, sd_int = input$sd_int, sd_slope = input$sd_slope,
    cor_is = input$cor_is, sd_resid = input$sd_resid,
    beta = c(input$betaA, input$betaB), lambda = input$lambda,
    depends_on = input$depends)))
  preview <- eventReactive(input$go_preview, {
    m <- appc_model()
    withProgress(message = "Calibrating the Appendix C model", value = 0.2, {
      ca <- safe(cs_calibrate_appc(m, "A", input$pA, input$muA, n_cal = input$n_cal))
      setProgress(0.6)
      cb <- safe(cs_calibrate_appc(m, "B", input$pB, input$muB, n_cal = input$n_cal))
      setProgress(0.8)
      pe <- cs_appc_policy_effect(m, ca, cb, n = input$n_cal)
    })
    list(m = m, ca = ca, cb = cb, pe = pe)
  })
  output$appc_plot <- renderPlot({
    p <- preview()
    cs_appc_plot(p$m, p$ca, p$cb, n = min(input$n_cal, 2e4))
  })
  output$appc_txt <- renderPrint({
    p <- preview()
    for (cal in list(p$ca, p$cb))
      cat(sprintf("Arm %s: alpha = %.3f, final mean on treatment = %.2f, completion %.3f, completers' mean %.2f (SD %.2f)\n",
                  cal$arm, cal$alpha, cal$final_mean, cal$completion,
                  cal$completers_mean, cal$completers_sd))
    cat(sprintf("Hypothetical (on-treatment) effect: %.3f; treatment policy effect: %.3f; completers' effect: %.3f\n",
                p$ca$final_mean - p$cb$final_mean, p$pe, input$muA - input$muB))
  })
  sim <- eventReactive(input$go_sim, {
    d <- cs_des()
    validate(need(input$estimand == "cs" && d$efficacy == "continuous" &&
                    d$margin_d == 0 && d$margin_e == 0,
                  "Simulation needs the composite test with a continuous outcome (superiority)."))
    n <- if (is.na(input$sim_n)) cs_sample_size(d, input$power)$nA else input$sim_n
    appc <- input$sim_dropout == "appendixC"
    withProgress(message = "Simulating trials", value = 0, {
      step <- max(1, input$sim_nsim %/% 50)
      safe(cs_simulate(d, n, nsim = input$sim_nsim, dropout = input$sim_dropout,
                       model = if (appc) appc_model() else NULL,
                       n_cal = input$n_cal, null = input$sim_null,
                       muB = input$muB, s = input$s,
                       comparators = if (appc) input$comparators else character(0),
                       M = input$M, mmrm_engine = input$engine,
                       seed = input$sim_seed,
                       progress = function(i, N) if (i %% step == 0)
                         setProgress(i / N, detail = sprintf("%d of %d", i, N))))
    })
  })
  output$sim_txt <- renderPrint(print(sim()))
  output$cmp_sim_tab <- renderTable({
    x <- sim()$comparators
    validate(need(!is.null(x), ""))
    x
  }, digits = 3)

  # ---- Analyze tab
  dat <- reactive({
    req(input$file)
    safe(utils::read.csv(input$file$datapath, stringsAsFactors = FALSE))
  })
  output$col_ui <- renderUI({
    nm <- names(dat())
    tagList(
      selectInput("c_arm", "Arm column", nm, nm[1]),
      uiOutput("active_ui"),
      selectInput("c_comp", "Completed (TRUE/FALSE or 1/0)", nm, nm[min(2, length(nm))]),
      selectInput("c_y", "Improvement (completers)", nm, nm[min(3, length(nm))]),
      selectInput("c_cov", "Baseline covariates", nm, multiple = TRUE),
      selectInput("c_strata", "Stratum", c("(none)", nm)))
  })
  output$active_ui <- renderUI({
    req(input$c_arm)
    lv <- sort(unique(dat()[[input$c_arm]]))
    selectInput("c_active", "Active arm", lv, lv[1])
  })
  output$an_head <- renderTable(utils::head(dat(), 8))
  an <- eventReactive(input$go_an, {
    req(input$c_active)
    wd <- if (input$an_method == "fisher") 0.5 else input$an_wd
    safe(cs_analyze(dat(), arm = input$c_arm, active = input$c_active,
                    completed = input$c_comp, outcome = input$c_y,
                    covariates = if (length(input$c_cov)) input$c_cov else NULL,
                    strata = if (input$c_strata == "(none)") NULL else input$c_strata,
                    weights = if (input$an_method == "wz" && isTRUE(input$an_aligned))
                      "aligned" else c(wd, 1 - wd),
                    method = input$an_method,
                    contingent = input$an_cont, c = input$an_c,
                    alpha = input$an_alpha, s = input$an_s,
                    psi_type = input$an_psi, bootstrap = input$an_boot, seed = 1))
  })
  output$an_txt <- renderPrint(print(an()))
}

shinyApp(ui, server)

