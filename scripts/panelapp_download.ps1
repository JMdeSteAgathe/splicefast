# Telecharge toutes les pages de l'API PanelApp Australie (genes) en JSON.
# Utilise le proxy et les certificats de Windows.
param([string]$OutDir = ".\panelapp_json",
      [string]$Api = "https://panelapp-aus.org/api/v1/genes/")
$ErrorActionPreference = "Stop"
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$url = $Api; $n = 0
while ($url) {
  $n++
  for ($try = 1; $try -le 5; $try++) {
    try { $page = Invoke-RestMethod -Uri $url -Headers @{accept = "application/json"} -UseDefaultCredentials; break }
    catch { if ($try -eq 5) { throw }; Start-Sleep -Seconds (2 * $try) }
  }
  $page | ConvertTo-Json -Depth 20 -Compress | Set-Content -Encoding UTF8 -Path (Join-Path $OutDir "page_$n.json")
  if ($n % 20 -eq 0) { Write-Host "  page $n" }
  $url = $page.next
}
Write-Host "PanelApp : $n pages dans $OutDir"
