class MaintenanceController < ApplicationController

  def index
    if !current_user.present?
      redirect_to root_path
    end
    @labmethod = Labmethod.new
    @profile = Profile.new
    @labmethods = Labmethod.all
    # binding.pry
  end

  def show

  end
end
