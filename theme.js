// ✅ Initialize
showStartScreen();

  function setTheme(mode) {
  document.documentElement.setAttribute("data-theme", mode);
  localStorage.setItem("theme", mode);
  updateLogo();
}

function toggleTheme() {
  const cur = document.documentElement.getAttribute("data-theme");
  setTheme(cur === "dark" ? "light" : "dark");
}

function loadTheme() {
  const saved = localStorage.getItem("theme") || (window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
  setTheme(saved);
}

function updateLogo() {
  const logo = document.querySelector(".logo");
  if (!logo) return;
  const mode = document.documentElement.getAttribute("data-theme");
  logo.style.opacity = "0";
  setTimeout(() => {
    logo.src = mode === "dark" ? "./logo1.png" : "./logo.png";
    logo.style.opacity = "1";
  }, 150);
}

window.addEventListener('load', loadTheme);


