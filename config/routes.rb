Rails.application.routes.draw do
  devise_for :users
  # For details on the DSL available within this file, see http://guides.rubyonrails.org/routing.html
  root("static_pages#index")

  namespace :api do
    namespace :v1 do
      resources :interests, only: [:index, :show]
      resources :profiles, only: [:index]
	  resources :newstuffs, only: [:index]
      resources :labmethods, only: [:index, :show]
    end
  end

  resources :maintenance, only: [:index, :show]
  resources :labmethods, only: [:new, :create, :destroy]
  resources :profiles, only: [:new, :create, :edit, :update, :destroy]
  resources :newstuffs, only: [:new, :create, :edit, :update, :destroy]
  resources :pastmembers, only: [:new, :create, :edit, :update, :destroy]
  get '*path', to: 'static_pages#index'
end
