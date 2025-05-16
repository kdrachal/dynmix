\name{oil}
\alias{oil}

\title{Crude Oil Data.}

\description{Selected data from oil market.}

\usage{data(oil)}

\format{\code{oil} is \code{\link[base]{matrix}} object such that columnwise are
\itemize{
\item \code{crudeoil$p_oil} -- average spot price of crude oil (Brent, Dubai and WTI) in USD per barrel
\item \code{crudeoil$cons} -- U.S. product supplied of crude oil and petroleum products in thousand barrels 
\item \code{crudeoil$r} -- U.S. 3-month treasury bill secondary market rate in \%
\item \code{crudeoil$stocks} -- U.S. share prices index, 2015=100
\item \code{crudeoil$ex_rate} -- U.S. real effective exchange rate index (broad basket), 2020=100
}
}

\details{The data are in monthly frequency. They cover the period between Feb, 1998 and Oct, 2024. cons, r, stocks and ex_rate are lagged one period back.
}

\source{The data are provided by Bank for International Settlements, Board of Governors of the Federal Reserve System, OECD, U.S. Energy Information Administration and World Bank. 

\url{https://www.bis.org}

\url{https://www.eia.gov}

\url{https://www.federalreserve.gov}

\url{https://www.oecd.org/en.html} 

\url{https://www.worldbank.org/ext/en/home}
}

\examples{
data(oil)
}

\references{
Bank for International Settlements, 2025. Effective exchange rates, BIS WS_EER 1.0 (data set). \url{https://data.bis.org/topics/EER/BIS\%2CWS_EER\%2C1.0/M.R.B.US}

Board of Governors of the Federal Reserve System, 2025. Selected interest rates. \url{https://www.federalreserve.gov/releases/h15/}

OECD, 2025. Share prices. \url{https://www.oecd.org/en/data/indicators/share-prices.html}

U.S. Energy Information Administration, 2025. Petroleum /& other liquids. \url{https://www.eia.gov/petroleum/data.php}

World Bank, 2025. Commodity markets. \url{https://www.worldbank.org/en/research/commodity-markets}
}
