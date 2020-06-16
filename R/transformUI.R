#' User interface for flow cytometry transformation functions.
#'
#' This Shiny gadget provides a user interface for the log10, Logicle and asinh
#' transforms in the \code{\link[flowCore]{flowCore}} package.
#'
#' For each parameter, the user can choose between no, log10, Logicle or asinh
#' transformation and adjust the parameters for the Logicle and asinh functions.
#' When pressing the Set button, the current transformation for the selected
#' channel is added to the \code{\link[flowCore]{transformList}} that is
#' returned when the App is closed by pressing the Done button.
#'
#' For the logicle transformation, the transformation parameters are preset to
#' the values estimated by \code{\link[flowCore]{estimateLogicle}} for the
#' chosen channel. After manually adjusting the parameters, the parameters can
#' be reset to the estimated values by pressing the Auto button.
#'
#' Because of the longer calculation time, the parameters for the asinh
#' transformation are only estimated when the user presses the Auto button. You
#' can check the calculation progress in the console. Estimation is done via
#' \code{\link[flowVS]{estParamFlowVS}} from the \code{\link[flowVS]{flowVS}}
#' package.Without pressing the Auto button, fixed default values are set.
#'
#' To speed up calculations and data display, the app only works with a random
#' subset of the provided data. The size of that subset is defined by the
#' parameter \code{sampleSize} that is passed to the function.
#'
#' Data can be provided as a \code{\link[flowCore]{flowFrame}} or a
#' \code{\link[flowCore]{flowSet}}. If a \code{flowSet} is passed to the
#' function, all measurements in the \code{flowSet} are combined to a single
#' \code{flowFrame}, and the random subset is drawn from the combined data.
#'
#' Via the \code{transformList} argument, the user can pass a previously created
#' \code{transformList} to the function to visualize, modify and extend the
#' list. Note that currently all transformations not supported by transformUI
#' will be removed.
#'
#' Not that the density plot does not update live upon manipulating the sliders
#' but only after the user releases the slider.
#'
#' The app returns a \code{transformList} that can be used used with the
#' \code{\link[flowCore]{transform}} function in \code{flowCore}.
#'
#'
#'
#' @param dat either a \code{flowFrame} or \code{flowSet}.
#' @param sampleSize size of the subset used to estimate the parameters and
#'   display in the plots. Default is 10'000.
#' @param transformList a previously created \code{transformList}.
#'
#' @return transformation list with the channels for which a transformation has
#'   been set, NULL if no transformations were defined.
#' @export
#'
#' @author Stephan Benke (\email{s.benke@@cytometry.uzh.ch})
#'
transformUI <- function(dat, sampleSize = 1e4, transformList = NULL) {

  # downsampling for faster display
  if (class(dat) == "flowSet") {
    tmp <- do.call(rbind, flowCore::fsApply(dat, flowCore::exprs, simplify = F))
    ff <- flowCore::flowFrame(tmp[sample(nrow(tmp), sampleSize),])
  } else if (class(dat) == "flowFrame") {
    ff <- dat[sample(nrow(dat), sampleSize),]
  } else {
    stop("Input must be a flowFrame or flowSet.")
  }

  channels <- colnames(ff)

  # UI elements
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Transformation settings"),
    miniUI:: miniContentPanel(
      shiny::fluidRow(
        shiny::column(4,
                      shiny::tags$div(
                        style = "margin-bottom: 10px;",
                        shiny::uiOutput("channel"),
                        shiny::actionButton("ch.prev", "Previous"),
                        shiny::actionButton("ch.next", "Next"),
                        shiny::actionButton("ch.set", "Set")
                      ),
                      shiny::radioButtons("trans",
                                          shiny::h4("Transformation"),
                                          choices = c("None", "Log", "Logicle", "asinh"),
                                          selected = "None"),
                      shiny::uiOutput("sliders")
        ),
        shiny::column(8,
                      shiny::plotOutput("histo", height = "300px"),
                      shiny::uiOutput("description")
        )
      )
    )
  )

  # reactive part
  server <- function(input, output, session) {

    ## help text to  the transformations
    output$description <- shiny::renderUI({
      if (input$trans == "Log") {
        shiny::tags$div(
          style = "margin-left: 20px;",
          shiny::wellPanel(
            shiny::HTML("<b>Logarithmic</b> transformation to the base of 10.")
          )
        )
      } else if (input$trans == "Logicle") {
        shiny::tags$div(style = "margin-left: 20px;",
                        shiny::wellPanel(
                          shiny::h4("Logicle transformation"),
                          shiny::p("Automatically estimated by flowCore's estimateLogicle function."),
                          shiny::p("From the help file:"),
                          shiny::p("Logicle transformation creates a subset of biexponentialTransform hyperbolic sine transformation functions that provides several advantages over linear/log transformations for display of flow cytometry data."),
                          shiny::p("w: w is the linearization width in asymptotic decades. w should be > 0 and determines the slope of transformation at zero. w can be estimated using the equation w=(m-log10(t/abs(r)))/2, where r is the most negative value to be included in the display"),
                          shiny::p("t: Top of the scale data value, e.g, 10000 for common 4 decade data or 262144 for a 18 bit data range. t should be greater than zero"),
                          shiny::p("m: m is the full width of the transformed display in asymptotic decades. m should be greater than zero"),
                          shiny::p("a: Additional negative range to be included in the display in asymptotic decades. Positive values of the argument brings additional negative input values into the transformed display viewing area. Default value is zero corresponding to a Standard logicle function."),
                          shiny::p("Reference: Parks D.R., Roederer M., Moore W.A.(2006) A new \"logicle\" display method avoids deceptive effects of logarithmic scaling for low signals and compensated data. CytometryA, 96(6):541-51.")
                        )
        )

      } else if (input$trans == "asinh") {
        shiny::tags$div(style = "margin-left: 20px;",
                        shiny::wellPanel(
                          shiny::h4("arcsinh Transformation"),
                          shiny::p("According to the help file, the function definition is x <- asinh(a+b'*x)+c"),
                          shiny::p("Note that here we use b = 1/b' which is the more common parametrization of the cofactor."),
                          shiny::p("Pressing the Auto button will estimate b using the estParamFlowVS function from flowVS. The calculation takes several seconds.")
                        )
        )
      }
    })

    ## set transformation parameters
    if (is.null(transformList)) { # empty list if no transformList passed to transformUI
      transParams <- shiny::reactiveVal(list())
    } else {
      transParams <- shiny::reactiveVal(extractTransformListParams(ff, transformList))
    }

    shiny::observeEvent(input$ch.set, {
      ls <- transParams()
      if (shiny::req(input$trans) %in% c("None", "Log")) {
        ls[[shiny::req(input$channel)]] <- list(trans = input$trans)
      } else if (shiny::req(input$trans) == "Logicle") {
        ls[[shiny::req(input$channel)]] <- list(trans = input$trans,
                                                w = input$logicle.w,
                                                t = input$logicle.t,
                                                m = input$logicle.m,
                                                a = input$logicle.a)
      } else if (shiny::req(input$trans) == "asinh") {
        ls[[shiny::req(input$channel)]] <- list(trans = input$trans,
                                                a = input$asinh.a,
                                                b = input$asinh.b,
                                                c = input$asinh.c)
      }
      transParams(ls)
    })

    ## channel selection drop-down menu
    output$channel <- shiny::renderUI({
      shiny::selectInput("channel",
                         label = NULL,
                         choices = channels)
    })

    ## previous channel button
    shiny::observeEvent(input$ch.prev, {
      current <- which(input$channel == channels)
      if (current != 1) {
        shiny::updateSelectInput(session, "channel", selected = channels[current-1])
      }
    })

    ## next channel button
    shiny::observeEvent(input$ch.next, {
      current <- which(input$channel == channels)
      if (current != length(channels)) {
        shiny::updateSelectInput(session, "channel", selected = channels[current+1])
      }
    })

    ## conditional parameter input sliders
    output$sliders <- shiny::renderUI({

      if (input$trans == "Logicle") {
        if (!is.null(transParams()[[input$channel]])) {
          if (transParams()[[input$channel]]$trans == "Logicle") {
            params <- transParams()[[input$channel]][c("w", "t", "m", "a")]
          } else {
            params <- estimateLogicleParams(ff, input$channel)
          }
        } else {
          params <- estimateLogicleParams(ff, input$channel)
        }

        shiny::tags$div(
          shiny::sliderInput("logicle.w", "w", 0, 2, params$w, step = 0.01),
          shiny::sliderInput("logicle.t", "t", 1e4, 1e6, params$t, step = 1),
          shiny::sliderInput("logicle.m", "m", 1, 10, params$m, step = 0.1),
          shiny::sliderInput("logicle.a", "a", 0, 2, params$a, step = 0.01),
          shiny::actionButton("logicle.auto", "Auto")
        )
      } else if (input$trans == "asinh") {
        if (!is.null(transParams()[[input$channel]])) {
          if (transParams()[[input$channel]]$trans == "asinh") {
            params <- transParams()[[input$channel]][c("a", "b", "c")]
          } else {
            params <- list(a = 1, b = 9, c = 0)
          }
        } else {
          params <- list(a = 1, b = 9, c = 0)
        }

        shiny::tags$div(
          shiny::sliderInput("asinh.a", "a", 0, 2, params$a, step = 0.01),
          shiny::numericInput("asinh.b", "b", params$b, min = 0.1, step = 0.1),
          shiny::sliderInput("asinh.c", "c", 0, 10, params$c, step = 0.1),
          shiny::actionButton("asinh.auto", "Auto")
        )
      }
    })

    ## automatic logicle parameter estimation with estimateLogicle
    shiny::observeEvent(input$logicle.auto, {
      tryCatch( # to handle cases where estimateLogicle fails
        est <- flowCore::estimateLogicle(ff, input$channel),
        error = function(e) {
          warning(e)
          warning("Reverting to default parameters.")
        },
        finally = {
          if (!is.null(est)) {
            params <- list(w = unname(environment(est@transforms[[1]]@f)$w),
                           t = environment(est@transforms[[1]]@f)$t,
                           m = environment(est@transforms[[1]]@f)$m,
                           a = environment(est@transforms[[1]]@f)$a)
          } else {
            params <- list(w = 0.5, t = 262144, m = 4.5, a = 0)
          }
        }
      )
      shiny::updateSliderInput(session, "logicle.w", value = params$w)
      shiny::updateSliderInput(session, "logicle.t", value = params$t)
      shiny::updateSliderInput(session, "logicle.m", value = params$m)
      shiny::updateSliderInput(session, "logicle.a", value = params$a)
    })

    ## automatic asinh cofactor estimation with flowVS
    shiny::observeEvent(input$asinh.auto, {
      tryCatch( # to handle cases where estimation fails
        b <- flowVS::estParamFlowVS(flowCore::flowSet(ff), channels = input$channel),
        error = function(e) {
          warning(e)
          warning("Reverting to default parameters.")
        },
        finally = {
          if (is.null(b)) {
            b <- 9
          }
        }
      )
      shiny::updateNumericInput(session, "asinh.b", value = b)
    })

    ## transformation type selection based on channel type
    shiny::observeEvent(input$channel, {
      if (is.null(transParams()[[input$channel]])) {
        if (stringr::str_detect(input$channel, stringr::fixed("SC")) || stringr::str_detect(input$channel, stringr::fixed("Time"))) {
          shiny::updateRadioButtons(session, "trans", selected = "None")
        } else {
          shiny::updateRadioButtons(session, "trans", selected = "Logicle")
        }
      } else {
        shiny::updateRadioButtons(session, "trans", selected = transParams()[[input$channel]]$trans)
      }
    })

    ## data transformation
    dat <- shiny::reactive({
      if (input$trans == "None") {
        tibble::as_tibble(flowCore::exprs(ff))
      } else if (input$trans == "Log") {
        tibble::as_tibble(flowCore::exprs(flowCore::transform(ff, flowCore::transformList(shiny::req(input$channel),
                                                                                          flowCore::logTransform()))))
      } else if (input$trans == "Logicle") {
        tibble::as_tibble(flowCore::exprs(flowCore::transform(ff, flowCore::transformList(shiny::req(input$channel),
                                                                                          flowCore::logicleTransform(w = shiny::req(input$logicle.w),
                                                                                                                     t = shiny::req(input$logicle.t),
                                                                                                                     m = shiny::req(input$logicle.m),
                                                                                                                     a = shiny::req(input$logicle.a))))))
      } else if (input$trans == "asinh") {
        tibble::as_tibble(flowCore::exprs(flowCore::transform(ff, flowCore::transformList(shiny::req(input$channel),
                                                                                          flowCore::arcsinhTransform(a = shiny::req(input$asinh.a),
                                                                                                                     b = 1/shiny::req(input$asinh.b),
                                                                                                                     c = shiny::req(input$asinh.c))))))
      }
    })

    ## 1D density plot
    output$histo <- shiny::renderPlot({
      # scale_x_logicle seems to be broken atm, so we have to construct our own scale labels.

      breaks <- c(-1*10^seq(3, 2, -1), 0, 10^(2:7))
      if (input$trans == "None") {
        scale <- ggplot2::scale_x_continuous()
      } else if (input$trans == "Log") {
        scale <- ggplot2::scale_x_continuous(breaks = log10(breaks), labels = breaks)
      } else if (input$trans == "Logicle") {
        b <- flowCore::transform(flowCore::flowFrame(matrix(breaks, ncol = 1, dimnames = list(NULL, "a"))),
                                 flowCore::transformList("a", flowCore::logicleTransform(w = shiny::req(input$logicle.w),
                                                                                         t = shiny::req(input$logicle.t),
                                                                                         m = shiny::req(input$logicle.m),
                                                                                         a = shiny::req(input$logicle.a))))@exprs[,1]
        scale <- ggplot2::scale_x_continuous(breaks = b, labels = breaks)
      } else if (input$trans == "asinh") {
        b <- flowCore::transform(flowCore::flowFrame(matrix(breaks, ncol = 1, dimnames = list(NULL, "a"))),
                                 flowCore::transformList("a", flowCore::arcsinhTransform(a = shiny::req(input$asinh.a),
                                                                                         b = 1/shiny::req(input$asinh.b),
                                                                                         c = shiny::req(input$asinh.c))))@exprs[,1]
        scale <- ggplot2::scale_x_continuous(breaks = b, labels = breaks)
      }

      ggplot2::ggplot(dat(), ggplot2::aes_(ggplot2::sym(shiny::req(input$channel)))) +
        ggplot2::geom_density() +
        scale +
        ggplot2::theme_bw()
    })

    # close app and return values
    shiny::observeEvent(input$done, {

      if (length(transParams()) == 0) { # no transformation set
        output <- NULL
      } else {

        channels <- purrr::map_chr(seq_along(transParams()), function(i) {
          if (transParams()[[i]]$trans == "None") {
            NA
          } else {
            names(transParams())[i]
          }
        })
        channels <- channels[!is.na(channels)]

        if (length(channels) == 0) { # only None transformations set
          output <- NULL
        } else {

          # create transformation list

          transforms <- purrr::map(channels, function(x) {
            if (transParams()[[x]]$trans == "Log") {
              flowCore::logTransform()
            } else if (transParams()[[x]]$trans == "Logicle") {
              # params must be assigned first, otherwise expression is not evaluated before returning!?
              w <- transParams()[[x]]$w
              t <- transParams()[[x]]$t
              m <- transParams()[[x]]$m
              a <- transParams()[[x]]$a
              flowCore::logicleTransform(w = w,
                                         t = t,
                                         m = m,
                                         a = a)
            } else if (transParams()[[x]]$trans == "asinh") {
              # params must be assigned first, otherwise expression is not evaluated before returning!?
              a <- transParams()[[x]]$a
              b <- 1/transParams()[[x]]$b
              c <- transParams()[[x]]$c
              flowCore::arcsinhTransform(a = a, b = b, c = c)
            }
          })

          output <- flowCore::transformList(channels, transforms)
        }
      }

      shiny::stopApp(output)
    })
  }

  # gadget start
  shiny::runGadget(ui, server, viewer = shiny::paneViewer())
}


#' Estimate parameters of Logicle Transform
#'
#' Estimates the transformation parameters using estimateLogicle and returns the
#' found parameters instead of a transformation list. Returns default parameters
#' and a warning if estimateLogicle fails.
#'
#' Intended to be used internally by transformUI.
#'
#' @param ff a flowFrame
#' @param channel string indicating the channel for which to estimate the
#'   transformation. Intended to run only for one channel at a time.
#'
#' @return list with parameter values w, t, m and a.
#'
estimateLogicleParams <- function(ff, channel) {
  tryCatch( # to handle cases where estimateLogicle fails
    est <- flowCore::estimateLogicle(ff, channel),
    error = function(e) {
      warning(e)
      warning("Reverting to default parameters.")
    },
    finally = {
      if (!is.null(est)) {
        params <- list(w = unname(environment(est@transforms[[1]]@f)$w),
                       t = environment(est@transforms[[1]]@f)$t,
                       m = environment(est@transforms[[1]]@f)$m,
                       a = environment(est@transforms[[1]]@f)$a)
      } else {
        params <- list(w = 0.5, t = 262144, m = 4.5, a = 0)
      }
    }
  )
  params
}


#' Extract transformation parameters from transform list
#'
#' Intended to be used internally by transformUI.
#'
#' @param ff flowFrame
#' @param tl transformList
#'
#' @return list with transformations and parameters per channel
extractTransformListParams <- function(ff, tl) {
  out <- list()

  for (i in seq_along(tl@transforms)) {

    if (tl@transforms[[i]]@input %in% colnames(ff)) {
      if (environment(tl@transforms[[i]]@f)$transformationId %in%
          c("defaultLogTransform", "defaultArcsinhTransform", "defaultLogicleTransform")) {
        if (environment(tl@transforms[[i]]@f)$transformationId == "defaultLogTransform") {
          out[[tl@transforms[[i]]@input]] <- list(trans = "Log")
        } else if (environment(tl@transforms[[i]]@f)$transformationId == "defaultArcsinhTransform") {
          out[[tl@transforms[[i]]@input]] <- list(trans = "asinh",
                                                  a = environment(tl@transforms[[i]]@f)$a,
                                                  b = 1/environment(tl@transforms[[i]]@f)$b,
                                                  c = environment(tl@transforms[[i]]@f)$c)
        } else if (environment(tl@transforms[[i]]@f)$transformationId == "defaultLogicleTransform") {
          out[[tl@transforms[[i]]@input]] <- list(trans = "Logicle",
                                                  w = environment(tl@transforms[[i]]@f)$w,
                                                  t = environment(tl@transforms[[i]]@f)$t,
                                                  m = environment(tl@transforms[[i]]@f)$m,
                                                  a = environment(tl@transforms[[i]]@f)$a)
        }
      } else {
        warning(paste0("Channel ", tl@transforms[[i]]@input, ": ",
                       environment(tl@transforms[[i]]@f)$transformationId,
                       " not yet supported by transformUI, reverting to None."))
      }
    } else {
      warning(paste0("Channel ", tl@transforms[[i]]@input, " not found in flowFrame. Transformation ignored."))
    }
  }

  out
}

