# CLAUDE.md

## Проверка изменений

- Все проверки по возможности проводить с реальным запуском кода (тесты, примеры,
  скрипты из `verification/`, бенчмарки), а не только чтением исходников.
- Для запуска использовать последнюю стабильную версию Julia. Если Julia в окружении
  нет или версия устаревшая — установить актуальную (например, через `juliaup`:
  `curl -fsSL https://install.julialang.org | sh -s -- --yes`, затем
  `juliaup add release && juliaup default release`). В облачных сессиях это делает
  `.claude/hooks/session-start.sh`.

## Команды

- Быстрые тесты: `julia --project=. -e 'using Pkg; Pkg.test()'`.
  Отдельные файлы: `Pkg.test(test_args = ["golden", "contracts"])` — имена без
  `test_` и `.jl`, список в `test/runtests.jl`.
- Расширенная верификация (как в CI на PR):
  `VASILEK_EXTENDED=1 julia --project=. -e 'using Pkg; Pkg.test()'`.
- Исследования: один раз `julia --project=verification -e 'using Pkg; Pkg.instantiate()'`,
  затем `julia --project=verification verification/<study>.jl` (нужна Julia ≥ 1.11).
- Документация: `julia --project=docs docs/make.jl`. Примеры из README и docs
  проверяются тестами (`test_readme.jl`, `test_docs.jl`).

## Совместимость

- Нижняя граница — Julia 1.10 (LTS); CI гоняет `lts` и `1` на Linux, Windows и macOS.
  Работать на последней версии, но не использовать возможности языка и Pkg,
  которых нет в 1.10 (кроме окружений `verification/` и `docs/`, там пол — 1.11).
- Manifest не коммитится.

## Требования к коду

- Шаг солвера не выделяет память (`test/test_allocations.jl`) и работает в типе
  данных, а не только во `Float64`.
- Каждое заявление о точности из `verification/` закрепляется тестом; допуск и
  измеренное значение заносятся в таблицу `docs/src/verification.md`.
- Эталонные данные не правятся руками — перегенерируются через
  `test/generate_golden.jl`.
- Единицы и нормировка — `docs/src/normalization.md`.

## Оформление

- Каждый PR дописывает запись в `CHANGELOG.md` (Keep a Changelog, раздел текущей
  версии) с измерениями, которые обосновывают изменение. Ломающие изменения
  дополнительно описываются в `docs/src/migration-0.2.md`.
- Коммиты: `fix:` / `refactor:` и т. п., формулировка описывает поведение
  («a step of vlasov_poisson allocates nothing»).
- Ветку с master обновлять через merge, не через rebase.
- Код, комментарии, документация и коммиты — на английском; общение с
  пользователем — на русском.
