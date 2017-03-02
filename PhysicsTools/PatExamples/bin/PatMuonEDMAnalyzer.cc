





<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">



  <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/frameworks-eb63b5ca175f20151e478297e2205d5d6ac5acc8068a83c1d838e50e91689df4.css" media="all" rel="stylesheet" />
  <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/github-ce0c4ea874809df80d4e80fd270d0e13ea4f011c06ef9875b89fe1725aedff36.css" media="all" rel="stylesheet" />
  
  
  <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/site-f85b78c12cd707dafadc4d9e5e825e43bc3b3cc1283a814952649102189d669f.css" media="all" rel="stylesheet" />
  

  <meta name="viewport" content="width=device-width">
  
  <title>cmssw/PatMuonEDMAnalyzer.cc at tW · pphogat/cmssw · GitHub</title>
  <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub">
  <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub">
  <meta property="fb:app_id" content="1401488693436528">


  <link rel="assets" href="https://assets-cdn.github.com/">
  
  <meta name="pjax-timeout" content="1000">
  
  <meta name="request-id" content="B01A:5A8C:71C1CF:B7750D:58B7F6E4" data-pjax-transient>
  

  <meta name="selected-link" value="repo_source" data-pjax-transient>

  <meta name="google-site-verification" content="KT5gs8h0wvaagLKAVWq8bbeNwnZZK1r1XQysX3xurLU">
<meta name="google-site-verification" content="ZzhVyEFwb7w3e0-uOTltm8Jsck2F5StVihD0exw2fsA">
    <meta name="google-analytics" content="UA-3769691-2">

<meta content="collector.githubapp.com" name="octolytics-host" /><meta content="github" name="octolytics-app-id" /><meta content="https://collector.githubapp.com/github-external/browser_event" name="octolytics-event-url" /><meta content="B01A:5A8C:71C1CF:B7750D:58B7F6E4" name="octolytics-dimension-request_id" />
<meta content="/&lt;user-name&gt;/&lt;repo-name&gt;/blob/show" data-pjax-transient="true" name="analytics-location" />



  <meta class="js-ga-set" name="dimension1" content="Logged Out">



      <meta name="hostname" content="github.com">
  <meta name="user-login" content="">

      <meta name="expected-hostname" content="github.com">
    <meta name="js-proxy-site-detection-payload" content="Mzg4N2U0OTk5YmQ5NTAwODZmM2ZhM2Y5NjhjNjgyN2Y0OTYxM2NhODliN2Q4ZTAzZTU5ZTc0MGRkZWZiMDQ0MHx7InJlbW90ZV9hZGRyZXNzIjoiMTg4LjE4NC44OS4xNTMiLCJyZXF1ZXN0X2lkIjoiQjAxQTo1QThDOjcxQzFDRjpCNzc1MEQ6NThCN0Y2RTQiLCJ0aW1lc3RhbXAiOjE0ODg0NTEzMDAsImhvc3QiOiJnaXRodWIuY29tIn0=">


  <meta name="html-safe-nonce" content="476722059f7a48ba1c3bf53d70408053b324d5da">

  <meta http-equiv="x-pjax-version" content="e43befffed903825e51be774232c9c58">
  

    
  <meta name="description" content="cmssw - CMS Offline Software">
  <meta name="go-import" content="github.com/pphogat/cmssw git https://github.com/pphogat/cmssw.git">

  <meta content="12329789" name="octolytics-dimension-user_id" /><meta content="pphogat" name="octolytics-dimension-user_login" /><meta content="37315585" name="octolytics-dimension-repository_id" /><meta content="pphogat/cmssw" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="true" name="octolytics-dimension-repository_is_fork" /><meta content="10969551" name="octolytics-dimension-repository_parent_id" /><meta content="cms-sw/cmssw" name="octolytics-dimension-repository_parent_nwo" /><meta content="10969551" name="octolytics-dimension-repository_network_root_id" /><meta content="cms-sw/cmssw" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/pphogat/cmssw/commits/tW.atom" rel="alternate" title="Recent Commits to cmssw:tW" type="application/atom+xml">


    <link rel="canonical" href="https://github.com/pphogat/cmssw/blob/tW/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc" data-pjax-transient>


  <meta name="browser-stats-url" content="https://api.github.com/_private/browser/stats">

  <meta name="browser-errors-url" content="https://api.github.com/_private/browser/errors">

  <link rel="mask-icon" href="https://assets-cdn.github.com/pinned-octocat.svg" color="#000000">
  <link rel="icon" type="image/x-icon" href="https://assets-cdn.github.com/favicon.ico">

<meta name="theme-color" content="#1e2327">



  </head>

  <body class="logged-out env-production vis-public fork page-blob">
    

  <div class="position-relative js-header-wrapper ">
    <a href="#start-of-content" tabindex="1" class="accessibility-aid js-skip-to-content">Skip to content</a>
    <div id="js-pjax-loader-bar" class="pjax-loader-bar"><div class="progress"></div></div>

    
    
    



          <header class="site-header js-details-container Details" role="banner">
  <div class="container-responsive">
    <a class="header-logo-invertocat" href="https://github.com/" aria-label="Homepage" data-ga-click="(Logged out) Header, go to homepage, icon:logo-wordmark">
      <svg aria-hidden="true" class="octicon octicon-mark-github" height="32" version="1.1" viewBox="0 0 16 16" width="32"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"/></svg>
    </a>

    <button class="btn-link float-right site-header-toggle js-details-target" type="button" aria-label="Toggle navigation">
      <svg aria-hidden="true" class="octicon octicon-three-bars" height="24" version="1.1" viewBox="0 0 12 16" width="18"><path fill-rule="evenodd" d="M11.41 9H.59C0 9 0 8.59 0 8c0-.59 0-1 .59-1H11.4c.59 0 .59.41.59 1 0 .59 0 1-.59 1h.01zm0-4H.59C0 5 0 4.59 0 4c0-.59 0-1 .59-1H11.4c.59 0 .59.41.59 1 0 .59 0 1-.59 1h.01zM.59 11H11.4c.59 0 .59.41.59 1 0 .59 0 1-.59 1H.59C0 13 0 12.59 0 12c0-.59 0-1 .59-1z"/></svg>
    </button>

    <div class="site-header-menu">
      <nav class="site-header-nav">
          <a href="/features" class="js-selected-navigation-item nav-item" data-ga-click="Header, click, Nav menu - item:features" data-selected-links="/features /features">
            Features
</a>          <a href="/explore" class="js-selected-navigation-item nav-item" data-ga-click="Header, click, Nav menu - item:explore" data-selected-links="/explore /trending /trending/developers /integrations /integrations/feature/code /integrations/feature/collaborate /integrations/feature/ship /showcases /explore">
            Explore
</a>        <a href="/pricing" class="js-selected-navigation-item nav-item" data-ga-click="Header, click, Nav menu - item:pricing" data-selected-links="/pricing /pricing">
          Pricing
</a>      </nav>

      <div class="site-header-actions">
          <div class="header-search scoped-search site-scoped-search js-site-search" role="search">
  <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/pphogat/cmssw/search" class="js-site-search-form" data-scoped-search-url="/pphogat/cmssw/search" data-unscoped-search-url="/search" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
    <label class="form-control header-search-wrapper js-chromeless-input-container">
      <div class="header-search-scope">This repository</div>
      <input type="text"
        class="form-control header-search-input js-site-search-focus js-site-search-field is-clearable"
        data-hotkey="s"
        name="q"
        placeholder="Search"
        aria-label="Search this repository"
        data-unscoped-placeholder="Search GitHub"
        data-scoped-placeholder="Search"
        autocapitalize="off">
    </label>
</form></div>


          <a class="text-bold site-header-link" href="/login?return_to=%2Fpphogat%2Fcmssw%2Fblob%2FtW%2FPhysicsTools%2FPatExamples%2Fbin%2FPatMuonEDMAnalyzer.cc" data-ga-click="(Logged out) Header, clicked Sign in, text:sign-in">Sign in</a>
            <span class="text-gray">or</span>
            <a class="text-bold site-header-link" href="/join?source=header-repo" data-ga-click="(Logged out) Header, clicked Sign up, text:sign-up">Sign up</a>
      </div>
    </div>
  </div>
</header>


  </div>

  <div id="start-of-content" class="accessibility-aid"></div>

    <div id="js-flash-container">
</div>



  <div role="main">
      <div itemscope itemtype="http://schema.org/SoftwareSourceCode">
    <div id="js-repo-pjax-container" data-pjax-container>
      


<div class="pagehead repohead instapaper_ignore readability-menu experiment-repo-nav">
  <div class="container repohead-details-container">

    

<ul class="pagehead-actions">

  <li>
      <a href="/login?return_to=%2Fpphogat%2Fcmssw"
    class="btn btn-sm btn-with-count tooltipped tooltipped-n"
    aria-label="You must be signed in to watch a repository" rel="nofollow">
    <svg aria-hidden="true" class="octicon octicon-eye" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M8.06 2C3 2 0 8 0 8s3 6 8.06 6C13 14 16 8 16 8s-3-6-7.94-6zM8 12c-2.2 0-4-1.78-4-4 0-2.2 1.8-4 4-4 2.22 0 4 1.8 4 4 0 2.22-1.78 4-4 4zm2-4c0 1.11-.89 2-2 2-1.11 0-2-.89-2-2 0-1.11.89-2 2-2 1.11 0 2 .89 2 2z"/></svg>
    Watch
  </a>
  <a class="social-count" href="/pphogat/cmssw/watchers"
     aria-label="1 user is watching this repository">
    1
  </a>

  </li>

  <li>
      <a href="/login?return_to=%2Fpphogat%2Fcmssw"
    class="btn btn-sm btn-with-count tooltipped tooltipped-n"
    aria-label="You must be signed in to star a repository" rel="nofollow">
    <svg aria-hidden="true" class="octicon octicon-star" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M14 6l-4.9-.64L7 1 4.9 5.36 0 6l3.6 3.26L2.67 14 7 11.67 11.33 14l-.93-4.74z"/></svg>
    Star
  </a>

    <a class="social-count js-social-count" href="/pphogat/cmssw/stargazers"
      aria-label="0 users starred this repository">
      0
    </a>

  </li>

  <li>
      <a href="/login?return_to=%2Fpphogat%2Fcmssw"
        class="btn btn-sm btn-with-count tooltipped tooltipped-n"
        aria-label="You must be signed in to fork a repository" rel="nofollow">
        <svg aria-hidden="true" class="octicon octicon-repo-forked" height="16" version="1.1" viewBox="0 0 10 16" width="10"><path fill-rule="evenodd" d="M8 1a1.993 1.993 0 0 0-1 3.72V6L5 8 3 6V4.72A1.993 1.993 0 0 0 2 1a1.993 1.993 0 0 0-1 3.72V6.5l3 3v1.78A1.993 1.993 0 0 0 5 15a1.993 1.993 0 0 0 1-3.72V9.5l3-3V4.72A1.993 1.993 0 0 0 8 1zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3 10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3-10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"/></svg>
        Fork
      </a>

    <a href="/pphogat/cmssw/network" class="social-count"
       aria-label="2101 users forked this repository">
      2,101
    </a>
  </li>
</ul>

    <h1 class="public ">
  <svg aria-hidden="true" class="octicon octicon-repo-forked" height="16" version="1.1" viewBox="0 0 10 16" width="10"><path fill-rule="evenodd" d="M8 1a1.993 1.993 0 0 0-1 3.72V6L5 8 3 6V4.72A1.993 1.993 0 0 0 2 1a1.993 1.993 0 0 0-1 3.72V6.5l3 3v1.78A1.993 1.993 0 0 0 5 15a1.993 1.993 0 0 0 1-3.72V9.5l3-3V4.72A1.993 1.993 0 0 0 8 1zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3 10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3-10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"/></svg>
  <span class="author" itemprop="author"><a href="/pphogat" class="url fn" rel="author">pphogat</a></span><!--
--><span class="path-divider">/</span><!--
--><strong itemprop="name"><a href="/pphogat/cmssw" data-pjax="#js-repo-pjax-container">cmssw</a></strong>

    <span class="fork-flag">
      <span class="text">forked from <a href="/cms-sw/cmssw">cms-sw/cmssw</a></span>
    </span>
</h1>

  </div>
  <div class="container">
    
<nav class="reponav js-repo-nav js-sidenav-container-pjax"
     itemscope
     itemtype="http://schema.org/BreadcrumbList"
     role="navigation"
     data-pjax="#js-repo-pjax-container">

  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a href="/pphogat/cmssw/tree/tW" class="js-selected-navigation-item selected reponav-item" data-hotkey="g c" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches /pphogat/cmssw/tree/tW" itemprop="url">
      <svg aria-hidden="true" class="octicon octicon-code" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M9.5 3L8 4.5 11.5 8 8 11.5 9.5 13 14 8 9.5 3zm-5 0L0 8l4.5 5L6 11.5 2.5 8 6 4.5 4.5 3z"/></svg>
      <span itemprop="name">Code</span>
      <meta itemprop="position" content="1">
</a>  </span>


  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a href="/pphogat/cmssw/pulls" class="js-selected-navigation-item reponav-item" data-hotkey="g p" data-selected-links="repo_pulls /pphogat/cmssw/pulls" itemprop="url">
      <svg aria-hidden="true" class="octicon octicon-git-pull-request" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M11 11.28V5c-.03-.78-.34-1.47-.94-2.06C9.46 2.35 8.78 2.03 8 2H7V0L4 3l3 3V4h1c.27.02.48.11.69.31.21.2.3.42.31.69v6.28A1.993 1.993 0 0 0 10 15a1.993 1.993 0 0 0 1-3.72zm-1 2.92c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zM4 3c0-1.11-.89-2-2-2a1.993 1.993 0 0 0-1 3.72v6.56A1.993 1.993 0 0 0 2 15a1.993 1.993 0 0 0 1-3.72V4.72c.59-.34 1-.98 1-1.72zm-.8 10c0 .66-.55 1.2-1.2 1.2-.65 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"/></svg>
      <span itemprop="name">Pull requests</span>
      <span class="counter">0</span>
      <meta itemprop="position" content="3">
</a>  </span>

  <a href="/pphogat/cmssw/projects" class="js-selected-navigation-item reponav-item" data-selected-links="repo_projects new_repo_project repo_project /pphogat/cmssw/projects">
    <svg aria-hidden="true" class="octicon octicon-project" height="16" version="1.1" viewBox="0 0 15 16" width="15"><path fill-rule="evenodd" d="M10 12h3V2h-3v10zm-4-2h3V2H6v8zm-4 4h3V2H2v12zm-1 1h13V1H1v14zM14 0H1a1 1 0 0 0-1 1v14a1 1 0 0 0 1 1h13a1 1 0 0 0 1-1V1a1 1 0 0 0-1-1z"/></svg>
    Projects
    <span class="counter">0</span>
</a>


  <a href="/pphogat/cmssw/pulse" class="js-selected-navigation-item reponav-item" data-selected-links="pulse /pphogat/cmssw/pulse">
    <svg aria-hidden="true" class="octicon octicon-pulse" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M11.5 8L8.8 5.4 6.6 8.5 5.5 1.6 2.38 8H0v2h3.6l.9-1.8.9 5.4L9 8.5l1.6 1.5H14V8z"/></svg>
    Pulse
</a>
  <a href="/pphogat/cmssw/graphs" class="js-selected-navigation-item reponav-item" data-selected-links="repo_graphs repo_contributors /pphogat/cmssw/graphs">
    <svg aria-hidden="true" class="octicon octicon-graph" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M16 14v1H0V0h1v14h15zM5 13H3V8h2v5zm4 0H7V3h2v10zm4 0h-2V6h2v7z"/></svg>
    Graphs
</a>

</nav>

  </div>
</div>

<div class="container new-discussion-timeline experiment-repo-nav">
  <div class="repository-content">

    

<a href="/pphogat/cmssw/blob/e9fc925f7a62163e1a7283136e8c58792b8a27f0/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc" class="d-none js-permalink-shortcut" data-hotkey="y">Permalink</a>

<!-- blob contrib key: blob_contributors:v21:d7601172d37650b95c7a6dbaa28b2bae -->

<div class="file-navigation js-zeroclipboard-container">
  
<div class="select-menu branch-select-menu js-menu-container js-select-menu float-left">
  <button class="btn btn-sm select-menu-button js-menu-target css-truncate" data-hotkey="w"
    
    type="button" aria-label="Switch branches or tags" tabindex="0" aria-haspopup="true">
    <i>Branch:</i>
    <span class="js-select-button css-truncate-target">tW</span>
  </button>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax aria-hidden="true">

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <svg aria-label="Close" class="octicon octicon-x js-menu-close" height="16" role="img" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
        <span class="select-menu-title">Switch branches/tags</span>
      </div>

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Filter branches/tags" id="context-commitish-filter-field" class="form-control js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" data-filter-placeholder="Filter branches/tags" class="js-select-menu-tab" role="tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" data-filter-placeholder="Find a tag…" class="js-select-menu-tab" role="tab">Tags</a>
            </li>
          </ul>
        </div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches" role="menu">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_4_1_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_4_1_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_4_1_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_4_2_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_4_2_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_4_2_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_4_4_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_4_4_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_4_4_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_5_2_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_5_2_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_5_2_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_5_3_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_5_3_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_5_3_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_6_1_X_SLHC/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_6_1_X_SLHC"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_6_1_X_SLHC
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_6_2_SLHCDEV_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_6_2_SLHCDEV_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_6_2_SLHCDEV_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_6_2_X_SLHC/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_6_2_X_SLHC"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_6_2_X_SLHC
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_6_2_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_6_2_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_6_2_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_6_2_0_SLHC12_patch/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_6_2_0_SLHC12_patch"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_6_2_0_SLHC12_patch
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_6_2_0_SLHC15_fixTk/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_6_2_0_SLHC15_fixTk"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_6_2_0_SLHC15_fixTk
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_0_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_0_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_0_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_0_XROOTD_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_0_XROOTD_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_0_XROOTD_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_0_0_pre/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_0_0_pre"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_0_0_pre
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_1_CLANG_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_1_CLANG_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_1_CLANG_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_1_DEVEL_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_1_DEVEL_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_1_DEVEL_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_1_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_1_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_1_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_1_4_patchX/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_1_4_patchX"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_1_4_patchX
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_1_9_patch/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_1_9_patch"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_1_9_patch
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_2_BOOSTIO_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_2_BOOSTIO_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_2_BOOSTIO_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_2_CLANG_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_2_CLANG_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_2_CLANG_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_2_DEVEL_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_2_DEVEL_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_2_DEVEL_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_2_GEANT10_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_2_GEANT10_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_2_GEANT10_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_2_ROOT6_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_2_ROOT6_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_2_ROOT6_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_2_THREADED_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_2_THREADED_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_2_THREADED_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_2_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_2_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_2_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_2_0_patch/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_2_0_patch"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_2_0_patch
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_2_0_pre3_conddb/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_2_0_pre3_conddb"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_2_0_pre3_conddb
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_3_CLANG_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_3_CLANG_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_3_CLANG_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_3_DEVEL_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_3_DEVEL_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_3_DEVEL_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_3_GEANT10_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_3_GEANT10_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_3_GEANT10_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_3_ROOT6_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_3_ROOT6_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_3_ROOT6_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_3_THREADED_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_3_THREADED_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_3_THREADED_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_3_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_3_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_3_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_CLANG_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_CLANG_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_CLANG_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_DEVEL_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_DEVEL_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_DEVEL_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_GEANT10_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_GEANT10_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_GEANT10_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_ROOT5_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_ROOT5_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_ROOT5_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_ROOT6_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_ROOT6_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_ROOT6_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_THREADED_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_THREADED_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_THREADED_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_0_pre/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_0_pre"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_0_pre
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_1_patchX/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_1_patchX"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_1_patchX
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_4_2_patchX/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_4_2_patchX"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_4_2_patchX
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_5_ROOT5_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_5_ROOT5_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_5_ROOT5_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_5_ROOT64_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_5_ROOT64_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_5_ROOT64_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/CMSSW_7_5_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="CMSSW_7_5_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                CMSSW_7_5_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/externaldecay-update-on-top-off-5_3_14/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="externaldecay-update-on-top-off-5_3_14"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                externaldecay-update-on-top-off-5_3_14
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/imported-CVS-HEAD/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="imported-CVS-HEAD"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                imported-CVS-HEAD
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/ktf-patch-1/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="ktf-patch-1"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                ktf-patch-1
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/ktf-patch-2/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="ktf-patch-2"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                ktf-patch-2
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-6751-NominalCollision2015-default-CMSSW_7_4_X/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-6751-NominalCollision2015-default-CMSSW_7_4_X"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-6751-NominalCollision2015-default-CMSSW_7_4_X
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-6776-740p1_taus_switchConfDB/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-6776-740p1_taus_switchConfDB"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-6776-740p1_taus_switchConfDB
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-6942-FSQDQM_74/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-6942-FSQDQM_74"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-6942-FSQDQM_74
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-7108-TriggerResultFilter_update_74x/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-7108-TriggerResultFilter_update_74x"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-7108-TriggerResultFilter_update_74x
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-7113-fixFillDescriptionInESRawToDigi/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-7113-fixFillDescriptionInESRawToDigi"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-7113-fixFillDescriptionInESRawToDigi
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-7613-revert_to_conddb1/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-7613-revert_to_conddb1"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-7613-revert_to_conddb1
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-7710-hcal_infinte_loop/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-7710-hcal_infinte_loop"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-7710-hcal_infinte_loop
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-8034-ExoValDev_74X_PR7_2/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-8034-ExoValDev_74X_PR7_2"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-8034-ExoValDev_74X_PR7_2
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-8504-develop_Validation_PostProcessor/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-8504-develop_Validation_PostProcessor"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-8504-develop_Validation_PostProcessor
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-9009-CMSSW_7_3_5_patch2/autoPseudoParabolic/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-9009-CMSSW_7_3_5_patch2/autoPseudoParabolic"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-9009-CMSSW_7_3_5_patch2/autoPseudoParabolic
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-9080-p5-addons/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-9080-p5-addons"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-9080-p5-addons
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-9212-labels_v2/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-9212-labels_v2"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-9212-labels_v2
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-9416-beamspot_kcanrebinFix_744p1/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-9416-beamspot_kcanrebinFix_744p1"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-9416-beamspot_kcanrebinFix_744p1
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/revert-9423-revert-9080-p5-addons/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="revert-9423-revert-9080-p5-addons"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-9423-revert-9080-p5-addons
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open selected"
               href="/pphogat/cmssw/blob/tW/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="tW"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                tW
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/test-commit/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="test-commit"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                test-commit
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/pphogat/cmssw/blob/test-new-prs/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
               data-name="test-new-prs"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                test-new-prs
              </span>
            </a>
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/untagged-215a6cfabb6564af9a73/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="untagged-215a6cfabb6564af9a73"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="untagged-215a6cfabb6564af9a73">
                untagged-215a6cfabb6564af9a73
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_0_pre5/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_0_pre5"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_0_pre5">
                CMSSW_7_5_0_pre5
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_0_pre4/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_0_pre4"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_0_pre4">
                CMSSW_7_5_0_pre4
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_0_pre4_ROOT5/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_0_pre4_ROOT5"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_0_pre4_ROOT5">
                CMSSW_7_5_0_pre4_ROOT5
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_0_pre3/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_0_pre3"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_0_pre3">
                CMSSW_7_5_0_pre3
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_0_pre3_ROOT5/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_0_pre3_ROOT5"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_0_pre3_ROOT5">
                CMSSW_7_5_0_pre3_ROOT5
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_0_pre2/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_0_pre2"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_0_pre2">
                CMSSW_7_5_0_pre2
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_0_pre1/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_0_pre1"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_0_pre1">
                CMSSW_7_5_0_pre1
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_0_pre1_ROOT5/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_0_pre1_ROOT5"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_0_pre1_ROOT5">
                CMSSW_7_5_0_pre1_ROOT5
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-12-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-12-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-12-1100">
                CMSSW_7_5_X_2015-06-12-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-11-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-11-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-11-2300">
                CMSSW_7_5_X_2015-06-11-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-11-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-11-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-11-1100">
                CMSSW_7_5_X_2015-06-11-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-10-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-10-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-10-2300">
                CMSSW_7_5_X_2015-06-10-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-10-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-10-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-10-1100">
                CMSSW_7_5_X_2015-06-10-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-09-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-09-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-09-2300">
                CMSSW_7_5_X_2015-06-09-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-09-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-09-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-09-1100">
                CMSSW_7_5_X_2015-06-09-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-08-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-08-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-08-2300">
                CMSSW_7_5_X_2015-06-08-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-08-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-08-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-08-1100">
                CMSSW_7_5_X_2015-06-08-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-07-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-07-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-07-2300">
                CMSSW_7_5_X_2015-06-07-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-06-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-06-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-06-2300">
                CMSSW_7_5_X_2015-06-06-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-06-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-06-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-06-1100">
                CMSSW_7_5_X_2015-06-06-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-05-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-05-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-05-2300">
                CMSSW_7_5_X_2015-06-05-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-05-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-05-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-05-1100">
                CMSSW_7_5_X_2015-06-05-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-04-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-04-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-04-2300">
                CMSSW_7_5_X_2015-06-04-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-04-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-04-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-04-1100">
                CMSSW_7_5_X_2015-06-04-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-03-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-03-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-03-2300">
                CMSSW_7_5_X_2015-06-03-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-03-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-03-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-03-1100">
                CMSSW_7_5_X_2015-06-03-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-02-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-02-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-02-2300">
                CMSSW_7_5_X_2015-06-02-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-02-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-02-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-02-1100">
                CMSSW_7_5_X_2015-06-02-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-01-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-01-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-01-2300">
                CMSSW_7_5_X_2015-06-01-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-06-01-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-06-01-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-06-01-1100">
                CMSSW_7_5_X_2015-06-01-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-31-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-31-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-31-2300">
                CMSSW_7_5_X_2015-05-31-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-31-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-31-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-31-1100">
                CMSSW_7_5_X_2015-05-31-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-30-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-30-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-30-2300">
                CMSSW_7_5_X_2015-05-30-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-30-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-30-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-30-1100">
                CMSSW_7_5_X_2015-05-30-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-29-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-29-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-29-2300">
                CMSSW_7_5_X_2015-05-29-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-29-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-29-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-29-1100">
                CMSSW_7_5_X_2015-05-29-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-28-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-28-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-28-2300">
                CMSSW_7_5_X_2015-05-28-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-28-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-28-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-28-1100">
                CMSSW_7_5_X_2015-05-28-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-27-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-27-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-27-2300">
                CMSSW_7_5_X_2015-05-27-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-27-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-27-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-27-1100">
                CMSSW_7_5_X_2015-05-27-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-26-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-26-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-26-2300">
                CMSSW_7_5_X_2015-05-26-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-26-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-26-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-26-1100">
                CMSSW_7_5_X_2015-05-26-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-25-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-25-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-25-2300">
                CMSSW_7_5_X_2015-05-25-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-25-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-25-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-25-1100">
                CMSSW_7_5_X_2015-05-25-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-23-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-23-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-23-2300">
                CMSSW_7_5_X_2015-05-23-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-23-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-23-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-23-1100">
                CMSSW_7_5_X_2015-05-23-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-22-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-22-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-22-2300">
                CMSSW_7_5_X_2015-05-22-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-22-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-22-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-22-1100">
                CMSSW_7_5_X_2015-05-22-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-21-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-21-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-21-2300">
                CMSSW_7_5_X_2015-05-21-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-21-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-21-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-21-1100">
                CMSSW_7_5_X_2015-05-21-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-20-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-20-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-20-2300">
                CMSSW_7_5_X_2015-05-20-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-20-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-20-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-20-1100">
                CMSSW_7_5_X_2015-05-20-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-19-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-19-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-19-2300">
                CMSSW_7_5_X_2015-05-19-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-19-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-19-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-19-1100">
                CMSSW_7_5_X_2015-05-19-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-18-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-18-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-18-2300">
                CMSSW_7_5_X_2015-05-18-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-18-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-18-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-18-1100">
                CMSSW_7_5_X_2015-05-18-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-17-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-17-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-17-2300">
                CMSSW_7_5_X_2015-05-17-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-17-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-17-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-17-1100">
                CMSSW_7_5_X_2015-05-17-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-16-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-16-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-16-2300">
                CMSSW_7_5_X_2015-05-16-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-16-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-16-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-16-1100">
                CMSSW_7_5_X_2015-05-16-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-15-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-15-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-15-2300">
                CMSSW_7_5_X_2015-05-15-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-15-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-15-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-15-1100">
                CMSSW_7_5_X_2015-05-15-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-14-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-14-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-14-2300">
                CMSSW_7_5_X_2015-05-14-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-14-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-14-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-14-1100">
                CMSSW_7_5_X_2015-05-14-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-13-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-13-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-13-2300">
                CMSSW_7_5_X_2015-05-13-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-13-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-13-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-13-1100">
                CMSSW_7_5_X_2015-05-13-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-12-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-12-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-12-2300">
                CMSSW_7_5_X_2015-05-12-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-12-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-12-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-12-1100">
                CMSSW_7_5_X_2015-05-12-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-11-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-11-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-11-2300">
                CMSSW_7_5_X_2015-05-11-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-11-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-11-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-11-1100">
                CMSSW_7_5_X_2015-05-11-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-10-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-10-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-10-2300">
                CMSSW_7_5_X_2015-05-10-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-10-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-10-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-10-1100">
                CMSSW_7_5_X_2015-05-10-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-09-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-09-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-09-2300">
                CMSSW_7_5_X_2015-05-09-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-09-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-09-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-09-1100">
                CMSSW_7_5_X_2015-05-09-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-08-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-08-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-08-2300">
                CMSSW_7_5_X_2015-05-08-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-08-1700/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-08-1700"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-08-1700">
                CMSSW_7_5_X_2015-05-08-1700
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-08-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-08-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-08-1100">
                CMSSW_7_5_X_2015-05-08-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-07-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-07-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-07-2300">
                CMSSW_7_5_X_2015-05-07-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-07-1200/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-07-1200"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-07-1200">
                CMSSW_7_5_X_2015-05-07-1200
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-07-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-07-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-07-1100">
                CMSSW_7_5_X_2015-05-07-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-06-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-06-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-06-2300">
                CMSSW_7_5_X_2015-05-06-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-06-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-06-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-06-1100">
                CMSSW_7_5_X_2015-05-06-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-05-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-05-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-05-2300">
                CMSSW_7_5_X_2015-05-05-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-05-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-05-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-05-1100">
                CMSSW_7_5_X_2015-05-05-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-04-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-04-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-04-2300">
                CMSSW_7_5_X_2015-05-04-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-04-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-04-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-04-1100">
                CMSSW_7_5_X_2015-05-04-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-03-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-03-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-03-2300">
                CMSSW_7_5_X_2015-05-03-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-03-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-03-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-03-1100">
                CMSSW_7_5_X_2015-05-03-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_X_2015-05-02-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_X_2015-05-02-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_X_2015-05-02-2300">
                CMSSW_7_5_X_2015-05-02-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-12-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-12-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-12-1100">
                CMSSW_7_5_THREADED_X_2015-06-12-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-11-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-11-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-11-2300">
                CMSSW_7_5_THREADED_X_2015-06-11-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-11-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-11-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-11-1100">
                CMSSW_7_5_THREADED_X_2015-06-11-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-10-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-10-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-10-2300">
                CMSSW_7_5_THREADED_X_2015-06-10-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-10-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-10-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-10-1100">
                CMSSW_7_5_THREADED_X_2015-06-10-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-09-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-09-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-09-2300">
                CMSSW_7_5_THREADED_X_2015-06-09-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-09-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-09-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-09-1100">
                CMSSW_7_5_THREADED_X_2015-06-09-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-08-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-08-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-08-2300">
                CMSSW_7_5_THREADED_X_2015-06-08-2300
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-08-1100/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-08-1100"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-08-1100">
                CMSSW_7_5_THREADED_X_2015-06-08-1100
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/pphogat/cmssw/tree/CMSSW_7_5_THREADED_X_2015-06-07-2300/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc"
              data-name="CMSSW_7_5_THREADED_X_2015-06-07-2300"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="CMSSW_7_5_THREADED_X_2015-06-07-2300">
                CMSSW_7_5_THREADED_X_2015-06-07-2300
              </span>
            </a>
        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div>

    </div>
  </div>
</div>

  <div class="BtnGroup float-right">
    <a href="/pphogat/cmssw/find/tW"
          class="js-pjax-capture-input btn btn-sm BtnGroup-item"
          data-pjax
          data-hotkey="t">
      Find file
    </a>
    <button aria-label="Copy file path to clipboard" class="js-zeroclipboard btn btn-sm BtnGroup-item tooltipped tooltipped-s" data-copied-hint="Copied!" type="button">Copy path</button>
  </div>
  <div class="breadcrumb js-zeroclipboard-target">
    <span class="repo-root js-repo-root"><span class="js-path-segment"><a href="/pphogat/cmssw/tree/tW"><span>cmssw</span></a></span></span><span class="separator">/</span><span class="js-path-segment"><a href="/pphogat/cmssw/tree/tW/PhysicsTools"><span>PhysicsTools</span></a></span><span class="separator">/</span><span class="js-path-segment"><a href="/pphogat/cmssw/tree/tW/PhysicsTools/PatExamples"><span>PatExamples</span></a></span><span class="separator">/</span><span class="js-path-segment"><a href="/pphogat/cmssw/tree/tW/PhysicsTools/PatExamples/bin"><span>bin</span></a></span><span class="separator">/</span><strong class="final-path">PatMuonEDMAnalyzer.cc</strong>
  </div>
</div>


  <div class="commit-tease">
      <span class="float-right">
        <a class="commit-tease-sha" href="/pphogat/cmssw/commit/e9fc925f7a62163e1a7283136e8c58792b8a27f0" data-pjax>
          e9fc925
        </a>
        <relative-time datetime="2017-03-02T10:30:43Z">Mar 2, 2017</relative-time>
      </span>
      <div>
        <img alt="" class="avatar" data-canonical-src="https://2.gravatar.com/avatar/177fee0bf60a74c027485992f88c46e3?d=https%3A%2F%2Fassets-cdn.github.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png&amp;r=x&amp;s=140" height="20" src="https://camo.githubusercontent.com/99cd461d3ff5cf1ea4f08668378747323d03ceb5/68747470733a2f2f322e67726176617461722e636f6d2f6176617461722f31373766656530626636306137346330323734383539393266383863343665333f643d68747470732533412532462532466173736574732d63646e2e6769746875622e636f6d253246696d6167657325324667726176617461727325324667726176617461722d757365722d3432302e706e6726723d7826733d313430" width="20" />
        <span class="user-mention">priyanka priyanka</span>
          <a href="/pphogat/cmssw/commit/e9fc925f7a62163e1a7283136e8c58792b8a27f0" class="message" data-pjax="true" title="changes done">changes done</a>
      </div>

    <div class="commit-tease-contributors">
      <button type="button" class="btn-link muted-link contributors-toggle" data-facebox="#blob_contributors_box">
        <strong>2</strong>
         contributors
      </button>
          <a class="avatar-link tooltipped tooltipped-s" aria-label="wmtan" href="/pphogat/cmssw/commits/tW/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc?author=wmtan"><img alt="@wmtan" class="avatar" height="20" src="https://avatars0.githubusercontent.com/u/4268350?v=3&amp;s=40" width="20" /> </a>
    <a class="avatar-link tooltipped tooltipped-s" aria-label="ktf" href="/pphogat/cmssw/commits/tW/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc?author=ktf"><img alt="@ktf" class="avatar" height="20" src="https://avatars0.githubusercontent.com/u/10544?v=3&amp;s=40" width="20" /> </a>


    </div>

    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header" data-facebox-id="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list" data-facebox-id="facebox-description">
          <li class="facebox-user-list-item">
            <img alt="@wmtan" height="24" src="https://avatars2.githubusercontent.com/u/4268350?v=3&amp;s=48" width="24" />
            <a href="/wmtan">wmtan</a>
          </li>
          <li class="facebox-user-list-item">
            <img alt="@ktf" height="24" src="https://avatars2.githubusercontent.com/u/10544?v=3&amp;s=48" width="24" />
            <a href="/ktf">ktf</a>
          </li>
      </ul>
    </div>
  </div>


<div class="file">
  <div class="file-header">
  <div class="file-actions">

    <div class="BtnGroup">
      <a href="/pphogat/cmssw/raw/tW/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc" class="btn btn-sm BtnGroup-item" id="raw-url">Raw</a>
        <a href="/pphogat/cmssw/blame/tW/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc" class="btn btn-sm js-update-url-with-hash BtnGroup-item" data-hotkey="b">Blame</a>
      <a href="/pphogat/cmssw/commits/tW/PhysicsTools/PatExamples/bin/PatMuonEDMAnalyzer.cc" class="btn btn-sm BtnGroup-item" rel="nofollow">History</a>
    </div>


        <button type="button" class="btn-octicon disabled tooltipped tooltipped-nw"
          aria-label="You must be signed in to make or propose changes">
          <svg aria-hidden="true" class="octicon octicon-pencil" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M0 12v3h3l8-8-3-3-8 8zm3 2H1v-2h1v1h1v1zm10.3-9.3L12 6 9 3l1.3-1.3a.996.996 0 0 1 1.41 0l1.59 1.59c.39.39.39 1.02 0 1.41z"/></svg>
        </button>
        <button type="button" class="btn-octicon btn-octicon-danger disabled tooltipped tooltipped-nw"
          aria-label="You must be signed in to make or propose changes">
          <svg aria-hidden="true" class="octicon octicon-trashcan" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M11 2H9c0-.55-.45-1-1-1H5c-.55 0-1 .45-1 1H2c-.55 0-1 .45-1 1v1c0 .55.45 1 1 1v9c0 .55.45 1 1 1h7c.55 0 1-.45 1-1V5c.55 0 1-.45 1-1V3c0-.55-.45-1-1-1zm-1 12H3V5h1v8h1V5h1v8h1V5h1v8h1V5h1v9zm1-10H2V3h9v1z"/></svg>
        </button>
  </div>

  <div class="file-info">
      112 lines (99 sloc)
      <span class="file-info-divider"></span>
    3.95 KB
  </div>
</div>

  

  <div itemprop="text" class="blob-wrapper data type-c">
      <table class="highlight tab-size js-file-line-container" data-tab-size="8">
      <tr>
        <td id="L1" class="blob-num js-line-number" data-line-number="1"></td>
        <td id="LC1" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>memory<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L2" class="blob-num js-line-number" data-line-number="2"></td>
        <td id="LC2" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>string<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L3" class="blob-num js-line-number" data-line-number="3"></td>
        <td id="LC3" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>vector<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L4" class="blob-num js-line-number" data-line-number="4"></td>
        <td id="LC4" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>sstream<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L5" class="blob-num js-line-number" data-line-number="5"></td>
        <td id="LC5" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>fstream<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L6" class="blob-num js-line-number" data-line-number="6"></td>
        <td id="LC6" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>iostream<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L7" class="blob-num js-line-number" data-line-number="7"></td>
        <td id="LC7" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L8" class="blob-num js-line-number" data-line-number="8"></td>
        <td id="LC8" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>TH1F.h<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L9" class="blob-num js-line-number" data-line-number="9"></td>
        <td id="LC9" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>TROOT.h<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L10" class="blob-num js-line-number" data-line-number="10"></td>
        <td id="LC10" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>TFile.h<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L11" class="blob-num js-line-number" data-line-number="11"></td>
        <td id="LC11" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>TSystem.h<span class="pl-pds">&gt;</span></span></td>
      </tr>
      <tr>
        <td id="L12" class="blob-num js-line-number" data-line-number="12"></td>
        <td id="LC12" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L13" class="blob-num js-line-number" data-line-number="13"></td>
        <td id="LC13" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>DataFormats/FWLite/interface/Event.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L14" class="blob-num js-line-number" data-line-number="14"></td>
        <td id="LC14" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>DataFormats/Common/interface/Handle.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L15" class="blob-num js-line-number" data-line-number="15"></td>
        <td id="LC15" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>FWCore/FWLite/interface/FWLiteEnabler.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L16" class="blob-num js-line-number" data-line-number="16"></td>
        <td id="LC16" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L17" class="blob-num js-line-number" data-line-number="17"></td>
        <td id="LC17" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>DataFormats/MuonReco/interface/Muon.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L18" class="blob-num js-line-number" data-line-number="18"></td>
        <td id="LC18" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>DataFormats/PatCandidates/interface/Muon.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L19" class="blob-num js-line-number" data-line-number="19"></td>
        <td id="LC19" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>PhysicsTools/FWLite/interface/TFileService.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L20" class="blob-num js-line-number" data-line-number="20"></td>
        <td id="LC20" class="blob-code blob-code-inner js-file-line">#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&quot;</span>PhysicsTools/FWLite/interface/CommandLineParser.h<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L21" class="blob-num js-line-number" data-line-number="21"></td>
        <td id="LC21" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L22" class="blob-num js-line-number" data-line-number="22"></td>
        <td id="LC22" class="blob-code blob-code-inner js-file-line"><span class="pl-k">int</span> <span class="pl-en">main</span>(<span class="pl-k">int</span> argc, <span class="pl-k">char</span>* argv[]) </td>
      </tr>
      <tr>
        <td id="L23" class="blob-num js-line-number" data-line-number="23"></td>
        <td id="LC23" class="blob-code blob-code-inner js-file-line">{</td>
      </tr>
      <tr>
        <td id="L24" class="blob-num js-line-number" data-line-number="24"></td>
        <td id="LC24" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> ----------------------------------------------------------------------</span></td>
      </tr>
      <tr>
        <td id="L25" class="blob-num js-line-number" data-line-number="25"></td>
        <td id="LC25" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> First Part: </span></td>
      </tr>
      <tr>
        <td id="L26" class="blob-num js-line-number" data-line-number="26"></td>
        <td id="LC26" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span></span></td>
      </tr>
      <tr>
        <td id="L27" class="blob-num js-line-number" data-line-number="27"></td>
        <td id="LC27" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span>  * enable FWLite </span></td>
      </tr>
      <tr>
        <td id="L28" class="blob-num js-line-number" data-line-number="28"></td>
        <td id="LC28" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span>  * book the histograms of interest </span></td>
      </tr>
      <tr>
        <td id="L29" class="blob-num js-line-number" data-line-number="29"></td>
        <td id="LC29" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span>  * open the input file</span></td>
      </tr>
      <tr>
        <td id="L30" class="blob-num js-line-number" data-line-number="30"></td>
        <td id="LC30" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> ----------------------------------------------------------------------</span></td>
      </tr>
      <tr>
        <td id="L31" class="blob-num js-line-number" data-line-number="31"></td>
        <td id="LC31" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L32" class="blob-num js-line-number" data-line-number="32"></td>
        <td id="LC32" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> load framework libraries</span></td>
      </tr>
      <tr>
        <td id="L33" class="blob-num js-line-number" data-line-number="33"></td>
        <td id="LC33" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">gSystem</span>-&gt;<span class="pl-c1">Load</span>( <span class="pl-s"><span class="pl-pds">&quot;</span>libFWCoreFWLite<span class="pl-pds">&quot;</span></span> );</td>
      </tr>
      <tr>
        <td id="L34" class="blob-num js-line-number" data-line-number="34"></td>
        <td id="LC34" class="blob-code blob-code-inner js-file-line">  <span class="pl-c1">FWLiteEnabler::enable</span>();</td>
      </tr>
      <tr>
        <td id="L35" class="blob-num js-line-number" data-line-number="35"></td>
        <td id="LC35" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L36" class="blob-num js-line-number" data-line-number="36"></td>
        <td id="LC36" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> initialize command line parser</span></td>
      </tr>
      <tr>
        <td id="L37" class="blob-num js-line-number" data-line-number="37"></td>
        <td id="LC37" class="blob-code blob-code-inner js-file-line">  optutl::CommandLineParser <span class="pl-smi">parser</span> (<span class="pl-s"><span class="pl-pds">&quot;</span>Analyze FWLite Histograms<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L38" class="blob-num js-line-number" data-line-number="38"></td>
        <td id="LC38" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L39" class="blob-num js-line-number" data-line-number="39"></td>
        <td id="LC39" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> set defaults</span></td>
      </tr>
      <tr>
        <td id="L40" class="blob-num js-line-number" data-line-number="40"></td>
        <td id="LC40" class="blob-code blob-code-inner js-file-line">  parser.<span class="pl-c1">integerValue</span> (<span class="pl-s"><span class="pl-pds">&quot;</span>maxEvents<span class="pl-pds">&quot;</span></span>  ) = <span class="pl-c1">1000</span>;</td>
      </tr>
      <tr>
        <td id="L41" class="blob-num js-line-number" data-line-number="41"></td>
        <td id="LC41" class="blob-code blob-code-inner js-file-line">  parser.<span class="pl-c1">integerValue</span> (<span class="pl-s"><span class="pl-pds">&quot;</span>outputEvery<span class="pl-pds">&quot;</span></span>) =   <span class="pl-c1">10</span>;</td>
      </tr>
      <tr>
        <td id="L42" class="blob-num js-line-number" data-line-number="42"></td>
        <td id="LC42" class="blob-code blob-code-inner js-file-line">  parser.<span class="pl-c1">stringValue</span>  (<span class="pl-s"><span class="pl-pds">&quot;</span>outputFile<span class="pl-pds">&quot;</span></span> ) = <span class="pl-s"><span class="pl-pds">&quot;</span>analyzeEdmTuple.root<span class="pl-pds">&quot;</span></span>;</td>
      </tr>
      <tr>
        <td id="L43" class="blob-num js-line-number" data-line-number="43"></td>
        <td id="LC43" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L44" class="blob-num js-line-number" data-line-number="44"></td>
        <td id="LC44" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> parse arguments</span></td>
      </tr>
      <tr>
        <td id="L45" class="blob-num js-line-number" data-line-number="45"></td>
        <td id="LC45" class="blob-code blob-code-inner js-file-line">  parser.<span class="pl-c1">parseArguments</span> (argc, argv);</td>
      </tr>
      <tr>
        <td id="L46" class="blob-num js-line-number" data-line-number="46"></td>
        <td id="LC46" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">int</span> maxEvents_ = parser.<span class="pl-c1">integerValue</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>maxEvents<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L47" class="blob-num js-line-number" data-line-number="47"></td>
        <td id="LC47" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">unsigned</span> <span class="pl-k">int</span> outputEvery_ = parser.<span class="pl-c1">integerValue</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>outputEvery<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L48" class="blob-num js-line-number" data-line-number="48"></td>
        <td id="LC48" class="blob-code blob-code-inner js-file-line">  std::string outputFile_ = parser.<span class="pl-c1">stringValue</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>outputFile<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L49" class="blob-num js-line-number" data-line-number="49"></td>
        <td id="LC49" class="blob-code blob-code-inner js-file-line">  std::vector&lt;std::string&gt; inputFiles_ = parser.<span class="pl-c1">stringVector</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>inputFiles<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L50" class="blob-num js-line-number" data-line-number="50"></td>
        <td id="LC50" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L51" class="blob-num js-line-number" data-line-number="51"></td>
        <td id="LC51" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> book a set of histograms</span></td>
      </tr>
      <tr>
        <td id="L52" class="blob-num js-line-number" data-line-number="52"></td>
        <td id="LC52" class="blob-code blob-code-inner js-file-line">  fwlite::TFileService fs = <span class="pl-c1">fwlite::TFileService</span>(outputFile_.<span class="pl-c1">c_str</span>());</td>
      </tr>
      <tr>
        <td id="L53" class="blob-num js-line-number" data-line-number="53"></td>
        <td id="LC53" class="blob-code blob-code-inner js-file-line">  TFileDirectory dir = fs.<span class="pl-c1">mkdir</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>analyzePatMuon<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L54" class="blob-num js-line-number" data-line-number="54"></td>
        <td id="LC54" class="blob-code blob-code-inner js-file-line">  TH1F* muonPt_  = dir.<span class="pl-smi">make</span>&lt;TH1F&gt;(<span class="pl-s"><span class="pl-pds">&quot;</span>muonPt<span class="pl-pds">&quot;</span></span>  , <span class="pl-s"><span class="pl-pds">&quot;</span>pt<span class="pl-pds">&quot;</span></span>  ,   <span class="pl-c1">100</span>,   <span class="pl-c1">0</span>., <span class="pl-c1">300</span>.);</td>
      </tr>
      <tr>
        <td id="L55" class="blob-num js-line-number" data-line-number="55"></td>
        <td id="LC55" class="blob-code blob-code-inner js-file-line">  TH1F* muonEta_ = dir.<span class="pl-smi">make</span>&lt;TH1F&gt;(<span class="pl-s"><span class="pl-pds">&quot;</span>muonEta<span class="pl-pds">&quot;</span></span> , <span class="pl-s"><span class="pl-pds">&quot;</span>eta<span class="pl-pds">&quot;</span></span> ,   <span class="pl-c1">100</span>,  -<span class="pl-c1">3</span>.,   <span class="pl-c1">3</span>.);</td>
      </tr>
      <tr>
        <td id="L56" class="blob-num js-line-number" data-line-number="56"></td>
        <td id="LC56" class="blob-code blob-code-inner js-file-line">  TH1F* muonPhi_ = dir.<span class="pl-smi">make</span>&lt;TH1F&gt;(<span class="pl-s"><span class="pl-pds">&quot;</span>muonPhi<span class="pl-pds">&quot;</span></span> , <span class="pl-s"><span class="pl-pds">&quot;</span>phi<span class="pl-pds">&quot;</span></span> ,   <span class="pl-c1">100</span>,  -<span class="pl-c1">5</span>.,   <span class="pl-c1">5</span>.);  </td>
      </tr>
      <tr>
        <td id="L57" class="blob-num js-line-number" data-line-number="57"></td>
        <td id="LC57" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L58" class="blob-num js-line-number" data-line-number="58"></td>
        <td id="LC58" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> loop the events</span></td>
      </tr>
      <tr>
        <td id="L59" class="blob-num js-line-number" data-line-number="59"></td>
        <td id="LC59" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">int</span> ievt=<span class="pl-c1">0</span>;  </td>
      </tr>
      <tr>
        <td id="L60" class="blob-num js-line-number" data-line-number="60"></td>
        <td id="LC60" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span>(<span class="pl-k">unsigned</span> <span class="pl-k">int</span> iFile=<span class="pl-c1">0</span>; iFile&lt;inputFiles_.<span class="pl-c1">size</span>(); ++iFile){</td>
      </tr>
      <tr>
        <td id="L61" class="blob-num js-line-number" data-line-number="61"></td>
        <td id="LC61" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">//</span> open input file (can be located on castor)</span></td>
      </tr>
      <tr>
        <td id="L62" class="blob-num js-line-number" data-line-number="62"></td>
        <td id="LC62" class="blob-code blob-code-inner js-file-line">    TFile* inFile = <span class="pl-c1">TFile::Open</span>(inputFiles_[iFile].<span class="pl-c1">c_str</span>());</td>
      </tr>
      <tr>
        <td id="L63" class="blob-num js-line-number" data-line-number="63"></td>
        <td id="LC63" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span>( inFile ){</td>
      </tr>
      <tr>
        <td id="L64" class="blob-num js-line-number" data-line-number="64"></td>
        <td id="LC64" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">//</span> ----------------------------------------------------------------------</span></td>
      </tr>
      <tr>
        <td id="L65" class="blob-num js-line-number" data-line-number="65"></td>
        <td id="LC65" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">//</span> Second Part: </span></td>
      </tr>
      <tr>
        <td id="L66" class="blob-num js-line-number" data-line-number="66"></td>
        <td id="LC66" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">//</span></span></td>
      </tr>
      <tr>
        <td id="L67" class="blob-num js-line-number" data-line-number="67"></td>
        <td id="LC67" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">//</span>  * loop the events in the input file </span></td>
      </tr>
      <tr>
        <td id="L68" class="blob-num js-line-number" data-line-number="68"></td>
        <td id="LC68" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">//</span>  * receive the collections of interest via fwlite::Handle</span></td>
      </tr>
      <tr>
        <td id="L69" class="blob-num js-line-number" data-line-number="69"></td>
        <td id="LC69" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">//</span>  * fill the histograms</span></td>
      </tr>
      <tr>
        <td id="L70" class="blob-num js-line-number" data-line-number="70"></td>
        <td id="LC70" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">//</span>  * after the loop close the input file</span></td>
      </tr>
      <tr>
        <td id="L71" class="blob-num js-line-number" data-line-number="71"></td>
        <td id="LC71" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">//</span> ----------------------------------------------------------------------      </span></td>
      </tr>
      <tr>
        <td id="L72" class="blob-num js-line-number" data-line-number="72"></td>
        <td id="LC72" class="blob-code blob-code-inner js-file-line">      fwlite::Event <span class="pl-smi">ev</span>(inFile);</td>
      </tr>
      <tr>
        <td id="L73" class="blob-num js-line-number" data-line-number="73"></td>
        <td id="LC73" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">for</span>(ev.<span class="pl-c1">toBegin</span>(); !ev.<span class="pl-c1">atEnd</span>(); ++ev, ++ievt){</td>
      </tr>
      <tr>
        <td id="L74" class="blob-num js-line-number" data-line-number="74"></td>
        <td id="LC74" class="blob-code blob-code-inner js-file-line">	edm::EventBase <span class="pl-k">const</span> &amp; event = ev;</td>
      </tr>
      <tr>
        <td id="L75" class="blob-num js-line-number" data-line-number="75"></td>
        <td id="LC75" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">//</span> break loop if maximal number of events is reached </span></td>
      </tr>
      <tr>
        <td id="L76" class="blob-num js-line-number" data-line-number="76"></td>
        <td id="LC76" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">if</span>(maxEvents_&gt;<span class="pl-c1">0</span> ? ievt+<span class="pl-c1">1</span>&gt;maxEvents_ : <span class="pl-c1">false</span>) <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L77" class="blob-num js-line-number" data-line-number="77"></td>
        <td id="LC77" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">//</span> simple event counter</span></td>
      </tr>
      <tr>
        <td id="L78" class="blob-num js-line-number" data-line-number="78"></td>
        <td id="LC78" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">if</span>(outputEvery_!=<span class="pl-c1">0</span> ? (ievt&gt;<span class="pl-c1">0</span> &amp;&amp; ievt%outputEvery_==<span class="pl-c1">0</span>) : <span class="pl-c1">false</span>) </td>
      </tr>
      <tr>
        <td id="L79" class="blob-num js-line-number" data-line-number="79"></td>
        <td id="LC79" class="blob-code blob-code-inner js-file-line">	  std::cout &lt;&lt; <span class="pl-s"><span class="pl-pds">&quot;</span>  processing event: <span class="pl-pds">&quot;</span></span> &lt;&lt; ievt &lt;&lt; std::endl;</td>
      </tr>
      <tr>
        <td id="L80" class="blob-num js-line-number" data-line-number="80"></td>
        <td id="LC80" class="blob-code blob-code-inner js-file-line"><span class="pl-c"><span class="pl-c">//</span>/////sdgkjsdfhnsdvdkjsdnfmdsbcjnsdbcnbvsx</span></td>
      </tr>
      <tr>
        <td id="L81" class="blob-num js-line-number" data-line-number="81"></td>
        <td id="LC81" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">//</span> Handle to the muon pt</span></td>
      </tr>
      <tr>
        <td id="L82" class="blob-num js-line-number" data-line-number="82"></td>
        <td id="LC82" class="blob-code blob-code-inner js-file-line">	edm::<span class="pl-c1">Handle</span>&lt;std::vector&lt;<span class="pl-k">float</span>&gt; &gt; muonPt;</td>
      </tr>
      <tr>
        <td id="L83" class="blob-num js-line-number" data-line-number="83"></td>
        <td id="LC83" class="blob-code blob-code-inner js-file-line">	event.<span class="pl-c1">getByLabel</span>(<span class="pl-c1">std::string</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>patMuonAnalyzer:pt<span class="pl-pds">&quot;</span></span>), muonPt);</td>
      </tr>
      <tr>
        <td id="L84" class="blob-num js-line-number" data-line-number="84"></td>
        <td id="LC84" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">//</span> loop muon collection and fill histograms</span></td>
      </tr>
      <tr>
        <td id="L85" class="blob-num js-line-number" data-line-number="85"></td>
        <td id="LC85" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">for</span>(std::vector&lt;<span class="pl-k">float</span>&gt;::const_iterator mu1=muonPt-&gt;<span class="pl-c1">begin</span>(); mu1!=muonPt-&gt;<span class="pl-c1">end</span>(); ++mu1){</td>
      </tr>
      <tr>
        <td id="L86" class="blob-num js-line-number" data-line-number="86"></td>
        <td id="LC86" class="blob-code blob-code-inner js-file-line">	  muonPt_ -&gt;<span class="pl-c1">Fill</span>( *mu1 );</td>
      </tr>
      <tr>
        <td id="L87" class="blob-num js-line-number" data-line-number="87"></td>
        <td id="LC87" class="blob-code blob-code-inner js-file-line">	}</td>
      </tr>
      <tr>
        <td id="L88" class="blob-num js-line-number" data-line-number="88"></td>
        <td id="LC88" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">//</span> Handle to the muon eta</span></td>
      </tr>
      <tr>
        <td id="L89" class="blob-num js-line-number" data-line-number="89"></td>
        <td id="LC89" class="blob-code blob-code-inner js-file-line">	edm::<span class="pl-c1">Handle</span>&lt;std::vector&lt;<span class="pl-k">float</span>&gt; &gt; muonEta;</td>
      </tr>
      <tr>
        <td id="L90" class="blob-num js-line-number" data-line-number="90"></td>
        <td id="LC90" class="blob-code blob-code-inner js-file-line">	event.<span class="pl-c1">getByLabel</span>(<span class="pl-c1">std::string</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>patMuonAnalyzer:eta<span class="pl-pds">&quot;</span></span>), muonEta);</td>
      </tr>
      <tr>
        <td id="L91" class="blob-num js-line-number" data-line-number="91"></td>
        <td id="LC91" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">for</span>(std::vector&lt;<span class="pl-k">float</span>&gt;::const_iterator mu1=muonEta-&gt;<span class="pl-c1">begin</span>(); mu1!=muonEta-&gt;<span class="pl-c1">end</span>(); ++mu1){</td>
      </tr>
      <tr>
        <td id="L92" class="blob-num js-line-number" data-line-number="92"></td>
        <td id="LC92" class="blob-code blob-code-inner js-file-line">	  muonEta_ -&gt;<span class="pl-c1">Fill</span>( *mu1 );</td>
      </tr>
      <tr>
        <td id="L93" class="blob-num js-line-number" data-line-number="93"></td>
        <td id="LC93" class="blob-code blob-code-inner js-file-line">	}</td>
      </tr>
      <tr>
        <td id="L94" class="blob-num js-line-number" data-line-number="94"></td>
        <td id="LC94" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">//</span> Handle to the muon phi</span></td>
      </tr>
      <tr>
        <td id="L95" class="blob-num js-line-number" data-line-number="95"></td>
        <td id="LC95" class="blob-code blob-code-inner js-file-line">	edm::<span class="pl-c1">Handle</span>&lt;std::vector&lt;<span class="pl-k">float</span>&gt; &gt; muonPhi;</td>
      </tr>
      <tr>
        <td id="L96" class="blob-num js-line-number" data-line-number="96"></td>
        <td id="LC96" class="blob-code blob-code-inner js-file-line">	event.<span class="pl-c1">getByLabel</span>(<span class="pl-c1">std::string</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>patMuonAnalyzer:phi<span class="pl-pds">&quot;</span></span>), muonPhi);</td>
      </tr>
      <tr>
        <td id="L97" class="blob-num js-line-number" data-line-number="97"></td>
        <td id="LC97" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">for</span>(std::vector&lt;<span class="pl-k">float</span>&gt;::const_iterator mu1=muonPhi-&gt;<span class="pl-c1">begin</span>(); mu1!=muonPhi-&gt;<span class="pl-c1">end</span>(); ++mu1){</td>
      </tr>
      <tr>
        <td id="L98" class="blob-num js-line-number" data-line-number="98"></td>
        <td id="LC98" class="blob-code blob-code-inner js-file-line">	  muonPhi_ -&gt;<span class="pl-c1">Fill</span>( *mu1 );</td>
      </tr>
      <tr>
        <td id="L99" class="blob-num js-line-number" data-line-number="99"></td>
        <td id="LC99" class="blob-code blob-code-inner js-file-line">	}</td>
      </tr>
      <tr>
        <td id="L100" class="blob-num js-line-number" data-line-number="100"></td>
        <td id="LC100" class="blob-code blob-code-inner js-file-line">      }  </td>
      </tr>
      <tr>
        <td id="L101" class="blob-num js-line-number" data-line-number="101"></td>
        <td id="LC101" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">//</span> close input file</span></td>
      </tr>
      <tr>
        <td id="L102" class="blob-num js-line-number" data-line-number="102"></td>
        <td id="LC102" class="blob-code blob-code-inner js-file-line">      inFile-&gt;<span class="pl-c1">Close</span>();</td>
      </tr>
      <tr>
        <td id="L103" class="blob-num js-line-number" data-line-number="103"></td>
        <td id="LC103" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L104" class="blob-num js-line-number" data-line-number="104"></td>
        <td id="LC104" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">//</span> break loop if maximal number of events is reached:</span></td>
      </tr>
      <tr>
        <td id="L105" class="blob-num js-line-number" data-line-number="105"></td>
        <td id="LC105" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">//</span> this has to be done twice to stop the file loop as well</span></td>
      </tr>
      <tr>
        <td id="L106" class="blob-num js-line-number" data-line-number="106"></td>
        <td id="LC106" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span>(maxEvents_&gt;<span class="pl-c1">0</span> ? ievt+<span class="pl-c1">1</span>&gt;maxEvents_ : <span class="pl-c1">false</span>) <span class="pl-k">break</span>;</td>
      </tr>
      <tr>
        <td id="L107" class="blob-num js-line-number" data-line-number="107"></td>
        <td id="LC107" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L108" class="blob-num js-line-number" data-line-number="108"></td>
        <td id="LC108" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">return</span> <span class="pl-c1">0</span>;</td>
      </tr>
      <tr>
        <td id="L109" class="blob-num js-line-number" data-line-number="109"></td>
        <td id="LC109" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L110" class="blob-num js-line-number" data-line-number="110"></td>
        <td id="LC110" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L111" class="blob-num js-line-number" data-line-number="111"></td>
        <td id="LC111" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
</table>

  </div>

</div>

<button type="button" data-facebox="#jump-to-line" data-facebox-class="linejump" data-hotkey="l" class="d-none">Jump to Line</button>
<div id="jump-to-line" style="display:none">
  <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="" class="js-jump-to-line-form" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
    <input class="form-control linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" aria-label="Jump to line" autofocus>
    <button type="submit" class="btn">Go</button>
</form></div>

  </div>
  <div class="modal-backdrop js-touch-events"></div>
</div>


    </div>
  </div>

  </div>

      <div class="container site-footer-container">
  <div class="site-footer" role="contentinfo">
    <ul class="site-footer-links float-right">
        <li><a href="https://github.com/contact" data-ga-click="Footer, go to contact, text:contact">Contact GitHub</a></li>
      <li><a href="https://developer.github.com" data-ga-click="Footer, go to api, text:api">API</a></li>
      <li><a href="https://training.github.com" data-ga-click="Footer, go to training, text:training">Training</a></li>
      <li><a href="https://shop.github.com" data-ga-click="Footer, go to shop, text:shop">Shop</a></li>
        <li><a href="https://github.com/blog" data-ga-click="Footer, go to blog, text:blog">Blog</a></li>
        <li><a href="https://github.com/about" data-ga-click="Footer, go to about, text:about">About</a></li>

    </ul>

    <a href="https://github.com" aria-label="Homepage" class="site-footer-mark" title="GitHub">
      <svg aria-hidden="true" class="octicon octicon-mark-github" height="24" version="1.1" viewBox="0 0 16 16" width="24"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"/></svg>
</a>
    <ul class="site-footer-links">
      <li>&copy; 2017 <span title="0.10780s from github-fe120-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="https://github.com/site/terms" data-ga-click="Footer, go to terms, text:terms">Terms</a></li>
        <li><a href="https://github.com/site/privacy" data-ga-click="Footer, go to privacy, text:privacy">Privacy</a></li>
        <li><a href="https://github.com/security" data-ga-click="Footer, go to security, text:security">Security</a></li>
        <li><a href="https://status.github.com/" data-ga-click="Footer, go to status, text:status">Status</a></li>
        <li><a href="https://help.github.com" data-ga-click="Footer, go to help, text:help">Help</a></li>
    </ul>
  </div>
</div>



  

  <div id="ajax-error-message" class="ajax-error-message flash flash-error">
    <svg aria-hidden="true" class="octicon octicon-alert" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M8.865 1.52c-.18-.31-.51-.5-.87-.5s-.69.19-.87.5L.275 13.5c-.18.31-.18.69 0 1 .19.31.52.5.87.5h13.7c.36 0 .69-.19.86-.5.17-.31.18-.69.01-1L8.865 1.52zM8.995 13h-2v-2h2v2zm0-3h-2V6h2v4z"/></svg>
    <button type="button" class="flash-close js-flash-close js-ajax-error-dismiss" aria-label="Dismiss error">
      <svg aria-hidden="true" class="octicon octicon-x" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
    </button>
    You can't perform that action at this time.
  </div>


    <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/compat-8e19569aacd39e737a14c8515582825f3c90d1794c0e5539f9b525b8eb8b5a8e.js"></script>
    <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/frameworks-506169cb2fe76254b921e8c944dd406e5cab2d719e657eace8ada98486231472.js"></script>
    <script async="async" crossorigin="anonymous" src="https://assets-cdn.github.com/assets/github-8b6edd55d17129e26ab9ec560b7e5ed2908e4babb81d02f8c229419dee58c0c1.js"></script>
    
    
    
    
  <div class="js-stale-session-flash stale-session-flash flash flash-warn flash-banner d-none">
    <svg aria-hidden="true" class="octicon octicon-alert" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M8.865 1.52c-.18-.31-.51-.5-.87-.5s-.69.19-.87.5L.275 13.5c-.18.31-.18.69 0 1 .19.31.52.5.87.5h13.7c.36 0 .69-.19.86-.5.17-.31.18-.69.01-1L8.865 1.52zM8.995 13h-2v-2h2v2zm0-3h-2V6h2v4z"/></svg>
    <span class="signed-in-tab-flash">You signed in with another tab or window. <a href="">Reload</a> to refresh your session.</span>
    <span class="signed-out-tab-flash">You signed out in another tab or window. <a href="">Reload</a> to refresh your session.</span>
  </div>
  <div class="facebox" id="facebox" style="display:none;">
  <div class="facebox-popup">
    <div class="facebox-content" role="dialog" aria-labelledby="facebox-header" aria-describedby="facebox-description">
    </div>
    <button type="button" class="facebox-close js-facebox-close" aria-label="Close modal">
      <svg aria-hidden="true" class="octicon octicon-x" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
    </button>
  </div>
</div>


  </body>
</html>

