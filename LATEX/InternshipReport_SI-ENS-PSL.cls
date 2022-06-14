% ---------------------------------------------------------------------
% Internship report M1S2
% ENS - March 2014
% Adapted from SelfArx class, 25/01/2012, v1.0
% ---------------------------------------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{InternshipReport-ENS-PSL}[29/03/2014, v1.0]
\RequirePackage{ifthen}
\RequirePackage{calc}
\AtEndOfClass{\RequirePackage{microtype}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*
\LoadClass{article}

\RequirePackage{times}      % Loads the Times-Roman Fonts
\RequirePackage{mathptmx}   % Loads the Times-Roman Math Fonts
\RequirePackage{ifpdf}      % Needed to pick between latex and pdflatex

% ---------------------------------------------------------------------
\RequirePackage[utf8]{inputenc}
\RequirePackage{amsmath,amsfonts,amssymb}
\RequirePackage{graphicx,xcolor}
\RequirePackage[english]{babel}
\RequirePackage{booktabs}
% ---------------------------------------------------------------------

\RequirePackage[version=3]{mhchem}

%%%%% for abstract+authors frames 

\graphicspath{{plots/}}

% ---------------------------------------------------------------------		  
% margins		  
\RequirePackage[left=2cm,%
                right=2cm,%
				top=2.25cm,%
				bottom=2.25cm,%
				headheight=11pt,%
				letterpaper]{geometry}%
\RequirePackage[labelfont={bf,sf},%
                labelsep=period,%
                justification=raggedright]{caption}
% ---------------------------------------------------------------------
\RequirePackage{fancyhdr}  % Needed to define custom headers/footers
\RequirePackage{lastpage}  % Number of pages in the document
\pagestyle{fancy}          % Enables the custom headers/footers
% Headers
\lhead{}%
\chead{}%
\rhead{\small\sffamily\bfseries\@ReportTitle\  --- \thepage/\pageref{LastPage}}
% Footers
\lfoot{}%
\cfoot{}%
\rfoot{}%
\renewcommand{\headrulewidth}{0pt}% % No header rule
\renewcommand{\footrulewidth}{0pt}% % No footer rule
% ---------------------------------------------------------------------
% section/subsection/paragraph set-up
\RequirePackage[explicit]{titlesec}
\titleformat{\section}
  {\color{color1}\large\sffamily\bfseries}
  {}
  {0em}
  {\colorbox{color2!10}{\parbox{\dimexpr\linewidth-2\fboxsep\relax}{\centering\arabic{section}. #1}}}
  []
\titleformat{name=\section,numberless}
  {\color{color1}\large\sffamily\bfseries}
  {}
  {0em}
  {\colorbox{color2!10}{\parbox{\dimexpr\linewidth-2\fboxsep\relax}{\centering#1}}}
  []  
\titleformat{\subsection}
  {\color{color1}\sffamily\bfseries}
  {\thesubsection}
  {0.5em}
  {#1}
  []
\titleformat{\subsubsection}
  {\sffamily\small\bfseries}
  {\thesubsubsection}
  {0.5em}
  {#1}
  []    
\titleformat{\paragraph}[runin]
  {\sffamily\small\bfseries}
  {}
  {0em}
  {#1} 
\titlespacing*{\section}{0pc}{3ex \@plus4pt \@minus3pt}{5pt}
\titlespacing*{\subsection}{0pc}{2.5ex \@plus3pt \@minus2pt}{0pt}
\titlespacing*{\subsubsection}{0pc}{2ex \@plus2.5pt \@minus1.5pt}{0pt}
\titlespacing*{\paragraph}{0pc}{1.5ex \@plus2pt \@minus1pt}{10pt}

% ---------------------------------------------------------------------
\RequirePackage{enumitem}
%\setlist{nolistsep} % Uncomment to remove spacing between items in lists (enumerate, itemize)
% ---------------------------------------------------------------------
% Remove brackets from numbering in List of References
\renewcommand{\@biblabel}[1]{\bfseries\color{color1}\textsuperscript{[#1]}}
%\setlength{\bibitemsep}{0cm}
% ---------------------------------------------------------------------
%\newcommand{\Year}[1]{\def\@Year{#1}}
\newcommand{\ReportTitle}[1]{\def\@ReportTitle{#1}}
\newcommand{\Author}[1]{\def\@Author{#1}}
\newcommand{\Supervisor}[1]{\def\@Supervisor{#1}}
\newcommand{\Laboratory}[1]{\def\@Laboratory{#1}}
\newcommand{\Abstract}[1]{\def\@Abstract{#1}}
\newcommand{\Keywords}[1]{\def\@Keywords{#1}}
% ---------------------------------------------------------------------
\renewcommand{\@maketitle}{%
\onecolumn%[{%
\thispagestyle{empty}%
\vskip-36pt%
\parbox{\textwidth-2\fboxsep-2\fboxrule}{%
\begin{minipage}[t][][b]{0.4\textwidth}\raggedright{\includegraphics[height=1.8cm]{plots/logo_ens_psl_couleur.png}}\end{minipage}
\hfill 
\begin{minipage}[t][][b]{0.4\textwidth}\raggedleft{
\includegraphics[trim={25cm 5cm 0 5cm},clip,height=1.3cm]{plots/ENS_Dpt_Chemistry_Logo.png}
\\ \vspace{0.5cm}
\sffamily\bfseries\selectfont Internship report M1 -- \the\year %\@Year
 }\end{minipage}\par}
\vskip20pt%
{\raggedright\color{color1}\sffamily\bfseries\fontsize{20}{25}\selectfont \@ReportTitle\par}%
\vskip10pt%
{\raggedright\color{color1}\sffamily\fontsize{12}{16}\selectfont  \@Author\par}%
\vskip10pt%
{\raggedright\sffamily\ {\bf \sffamily\textsc{supervisor}:} \@Supervisor\ - {\textit \@Laboratory\par}}%
\vskip18pt%
%\fcolorbox{color1}{white}{%
%\parbox{\textwidth-2\fboxsep-2\fboxrule}{\centering%
%\colorbox{color2!10}{%
%\parbox{\textwidth-3.5\fboxsep-3.5\fboxrule}{%
%\ifx\@Keywords\@empty
%\sffamily\small\textbf{\abstractname}\\\@Abstract
%\else
%\sffamily\small\textbf{\abstractname}\\\@Abstract\\[5pt]%
%\textbf{\keywordname}\\\@Keywords%
%\fi
%}%
%}%
%}%
%}%
\vskip25pt%
}%]%
}%
% ---------------------------------------------------------------------
\let\oldbibliography\thebibliography
\renewcommand{\thebibliography}[1]{%
\oldbibliography{#1}%
\setlength\itemsep{0pt}%
}
